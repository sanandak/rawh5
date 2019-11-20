package rawh5

import (
	"bytes"
	"compress/gzip"
	"database/sql"
	"encoding/binary"
	"encoding/hex"

	"fmt"
	"io"
	"io/ioutil"
	"os"
	"path/filepath"
	"strings"
	"time"

	_ "github.com/mattn/go-sqlite3"
	"github.com/sanandak/gp"
	"gonum.org/v1/hdf5"
)

const tformatStr = "20060102 150405" // go time string used for parsing

func RawH5(fname string, outdirname string) {

	pkt := make([]byte, 512)
	gpkt := gp.GPkt{}
	gps := make([]byte, 0)
	aux := make([]byte, 0)
	//ctr := make([]uint32, 0)
	//var avg [3]int
	var ns int
	fmt.Println("input file:", fname)
	fmt.Println("output directory", outdirname)

	// open input file
	file, err := os.Open(fname)
	if err != nil {
		fmt.Printf("file %s not opened err: %v\n", fname, err)
	}
	finfo, err := os.Stat(fname)
	if err != nil {
		fmt.Printf("file %s stat err: %v\n", fname, err)
	}
	defer file.Close()

	// npkts:number of pkts
	// nsmax:number of samps (this is an overestimate, b/c some pkts are gps/aux)
	npkts := finfo.Size() / 512
	nsmax := npkts * 30
	// metadata arrays
	gpsfix := make([]uint32, 0, npkts)
	ymd := make([]uint32, 0, npkts)
	hms := make([]uint32, 0, npkts)
	cpps := make([]uint32, 0, npkts)
	c0 := make([]uint32, 0, npkts)
	cN := make([]uint32, 0, npkts)
	tarr := make([]int64, 0, npkts)

	// array of all the samples interleaved (ch0,1,2,ctr)
	data := make([]int32, 0, nsmax*4)

	fmt.Println("npkts, nsmax", npkts, nsmax)

	// read pkts and parse
	var (
		t0                  time.Time
		oldT0               time.Time
		nseispkts           int64
		vernum              uint8
		sps                 uint
		dtsamp              float64 // samp int in s
		dtpktns             int64   // pkt interval + 1 samp in ns
		gainA, gainB, gainC float64
		euiStr              string
	)
	for {
		_, err = file.Read(pkt)
		if err == io.EOF {
			break
		}
		if gp.IsSeismic(pkt) {
			err = gp.Parsegp(pkt, &gpkt)
			if err != nil {
				continue
			}
			// 30 samps per pkt, 4 ch per samp...
			data = append(data, gpkt.Data[:]...)
			ns += 30

			// save metadata
			if nseispkts == 0 {
				//fmt.Printf("%+v\n", gpkt)
				vernum = gpkt.Vernum
				sps = gp.SPS[gpkt.SpsEnum]
				dtsamp = 1 / float64(sps)
				dtpktns = int64(dtsamp * 31 * 1e9)
				gA := (gpkt.GainABEnum >> 4) & 0x0f
				gB := gpkt.GainABEnum & 0x0f
				gC := (gpkt.GainCDEnum >> 4) & 0x0f
				//fmt.Println(gA, gB, gC)
				gainA = gp.Gain[gA]
				gainB = gp.Gain[gB]
				gainC = gp.Gain[gC]
				euiStr = strings.ToUpper(hex.EncodeToString([]byte{gpkt.Eui0, gpkt.Eui1, gpkt.Eui2}))
				fmt.Printf("ver %d sps %d dtsamp %f dtpktns %d gains %f %f %f eui %s\n",
					vernum, sps, dtsamp, dtpktns, gainA, gainB, gainC, euiStr)
			}

			tm, err := gp.GetSamp0T(gpkt)
			if err != nil {
				fmt.Println("err parsing gp...skipping", gpkt)
				continue
			}

			t0 = tm
			if dt := t0.Sub(oldT0).Nanoseconds(); dt > dtpktns || dt < -dtpktns {
				fmt.Printf("pkt %d time skip t0 %v oldt0 %v dt %d\n", nseispkts, t0, oldT0, dt)
			}
			oldT0 = t0
			tarr = append(tarr, t0.UnixNano())

			gpsfix = append(gpsfix, uint32(gpkt.Gpsfix))
			ymd = append(ymd, gpkt.YMD)
			hms = append(hms, gpkt.HMS)
			cpps = append(cpps, gpkt.Cpps)
			c0 = append(c0, gpkt.C0)
			cN = append(cN, gpkt.CN)

			nseispkts++
			if nseispkts%(npkts/10) == 0 {
				fmt.Printf("read %.0f%%\n", float64(nseispkts)/float64(npkts)*100)
				fmt.Printf("time %v %v\n", gpkt.YMD, gpkt.HMS)
			}

		} else if gp.IsGPS(pkt) {
			gps = append(gps, pkt...)
		} else if gp.IsAux(pkt) {
			aux = append(aux, pkt...)
		}
	}

	// get nodename
	nodename := euiStr //"UNKNOWN" //getNodename(euiStr, *dbname)
	fmt.Println("nodename = ", nodename)
	tStr := time.Unix(0, tarr[0]).UTC().Format("2006_01_02T15_04_05") // convert first time (a nanosec time) to time obj
	fmt.Println(tStr)
	// create hdf file and pebble group
	// fn: filename without directory or ext
	fn := fmt.Sprintf("%s_%s", nodename, tStr)
	gp5file := filepath.Join(outdirname, fn+".gp5")
	auxfile := filepath.Join(outdirname, fn+".aux")
	gpsfile := filepath.Join(outdirname, fn+".ubx")
	// write out ubx and aux buffers
	fmt.Println("saving to ubx file.gz", gpsfile)
	gzWrite(gps, gpsfile)
	fmt.Println("saving to aux file.gz", auxfile)
	gzWrite(aux, auxfile)
	fmt.Println("saving to gp5 file", gp5file)

	// create gp5file
	hfile, err := hdf5.CreateFile(gp5file, hdf5.F_ACC_TRUNC)
	if err != nil {
		panic("can't create file")
	}
	grp, err := hfile.CreateGroup(nodename)
	if err != nil {
		panic("can't create group")
	}
	defer grp.Close()
	defer hfile.Close()

	saveAttr("eui", grp, euiStr)
	saveAttr("sps", grp, float64(sps))
	saveAttr("gainA", grp, gainA)
	saveAttr("gainB", grp, gainB)
	saveAttr("gainC", grp, gainC)

	alldlen := uint(len(data))
	nslen := alldlen / 4
	fmt.Printf("dlen %d dlen/4 %d ns %d\n", alldlen, nslen, ns)

	// create data_a dataset (ns x 3)
	var fildims = []uint{alldlen / 4, 3}

	// chunking and compression
	var cdims = []uint{10000, 1}
	dcpl, err := hdf5.NewPropList(hdf5.P_DATASET_CREATE)
	if err != nil {
		panic(err)
	}
	err = dcpl.SetChunk(cdims)
	err = dcpl.SetDeflate(hdf5.DefaultCompression)

	filspace, err := hdf5.CreateSimpleDataspace(fildims, fildims)
	if err != nil {
		panic(err)
	}
	dset, err := grp.CreateDatasetWith("data_a", hdf5.T_NATIVE_INT, filspace, dcpl)
	if err != nil {
		panic(err)
	}
	defer dset.Close()
	filspace.Close()
	// save the three columns of data
	saveData(dset, data, 0)
	saveData(dset, data, 1)
	saveData(dset, data, 2)

	dcpl.Close()

	// save the metadata arrays
	saveMeta("ymd", grp, ymd)
	saveMeta("hms", grp, hms)
	saveMeta("gpsfix", grp, gpsfix)
	saveMeta("cpps", grp, cpps)
	saveMeta("c0", grp, c0)
	saveMeta("cN", grp, cN)
	saveTime("pktT0", grp, tarr)
	// // create dataset ctr_a (ns x 1)
	// // chunking and compression
	// cdims = []uint{10000}
	// dcpl, err = hdf5.NewPropList(hdf5.P_DATASET_CREATE)
	// if err != nil {
	// 	panic(err)
	// }
	// err = dcpl.SetChunk(cdims)
	// err = dcpl.SetDeflate(hdf5.DefaultCompression)

	// fildims = []uint{alldlen / 4}
	// filspace, err = hdf5.CreateSimpleDataspace(fildims, fildims)
	// if err != nil {
	// 	panic("can't create file dataspace for counter", err)
	// }
	// ctrdset, err := grp.CreateDatasetWith("ctr_a", hdf5.T_NATIVE_UINT, filspace, dcpl)
	// if err != nil {
	// 	panic("can't create dataset ctr_a", err)
	// }
	// defer ctrdset.Close()
	// filspace.Close()
	// saveCtr(ctrdset, data)
	dcpl.Close()

}

func gzWrite(b []byte, file string) {
	// compress data
	var zbuf bytes.Buffer
	var buf bytes.Buffer
	binary.Write(&buf, binary.LittleEndian, b)

	zw := gzip.NewWriter(&zbuf)
	if _, err := zw.Write(buf.Bytes()); err != nil {
		fmt.Println("err creating gzip buffer", err)
	}
	//fmt.Println("gz size", nwrt)
	zw.Close()

	// write out compressed data
	ioutil.WriteFile(file+".gz", zbuf.Bytes(), 0666)
}

func saveData(dset *hdf5.Dataset, data []int32, col uint) {
	// choose column of data = data_a(:, col)
	alldlen := uint(len(data))
	nslen := alldlen / 4
	offset := []uint{0, col}
	stride := []uint{1, 1}
	count := []uint{nslen, 1}
	block := []uint{1, 1}
	fildims := []uint{nslen, 3}
	memdims := []uint{alldlen}
	filspace, err := hdf5.CreateSimpleDataspace(fildims, fildims)
	if err != nil {
		panic(err)
		//panicf("can't create file dataspace for col %d %v", col, err)
	}
	defer filspace.Close()

	err = filspace.SelectHyperslab(offset, stride, count, block)
	if err != nil {
		panic(err)
		//panicf("can't create file dataspace hyperslab for col %d %v %v %v %v %v", col, offset, stride, count, block, err)
	}

	// choose the points in data corresponding to ch `col` - every fourth
	memspace, err := hdf5.CreateSimpleDataspace(memdims, memdims)
	if err != nil {
		panic(err)
		//panicf("can't create memory dataspace col %d, memdims: %v %v", col, memdims, err)
	}
	defer memspace.Close()

	moffset := []uint{0}
	mstride := []uint{4}
	mcount := []uint{nslen}
	mblock := []uint{1}
	memspace.SelectHyperslab(moffset, mstride, mcount, mblock)
	if err != nil {
		panic(err)
		//panicln("can't create memory hyperslab for col", col, moffset, mstride, mcount, mblock, err)
	}

	err = dset.WriteSubset(&data[0], memspace, filspace)
	if err != nil {
		panic(err)
		//panicln("can't write data_a col ", col, err)
	}

}

func saveCtr(dset *hdf5.Dataset, data []uint32) {
	// choose first/only column of data = ctr_a(:)

	alldlen := uint(len(data))
	nslen := alldlen / 4

	offset := []uint{0}
	stride := []uint{1}
	count := []uint{nslen}
	block := []uint{1}
	fildims := []uint{nslen}
	memdims := []uint{alldlen}
	filspace, err := hdf5.CreateSimpleDataspace(fildims, fildims)
	if err != nil {
		panic(err)
		//panic("can't create file dataspace", err)
	}
	defer filspace.Close()
	err = filspace.SelectHyperslab(offset, stride, count, block)
	if err != nil {
		panic(err)
		//panic("ctr: can't create file dataspace hyperslab", offset, stride, count, block, err)
	}
	// choose the points in data corresponding to ch 0 - every fourth
	memspace, err := hdf5.CreateSimpleDataspace(memdims, memdims)
	if err != nil {
		panic(err)
		//panic("can't create memory dataspace 0, memdims:", memdims, err)
	}
	defer memspace.Close()

	moffset := []uint{3}
	mstride := []uint{4}
	mcount := []uint{nslen}
	mblock := []uint{1}
	memspace.SelectHyperslab(moffset, mstride, mcount, mblock)
	if err != nil {
		panic(err)
		//panic("can't create memory hyperslab 0", moffset, mstride, mcount, mblock, err)
	}

	err = dset.WriteSubset(&data[0], memspace, filspace)
	if err != nil {
		panic(err)
		//panicln("can't write ctr_a", err)
	}
}

func saveMeta(name string, grp *hdf5.Group, arr []uint32) {
	nmeta := uint(len(arr))

	fildims := []uint{nmeta}
	filspace, err := hdf5.CreateSimpleDataspace(fildims, fildims)
	if err != nil {
		fmt.Println("can't create meta file dataspace", err)
	}
	defer filspace.Close()
	dset, err := grp.CreateDataset(name, hdf5.T_NATIVE_UINT, filspace)
	if err != nil {
		fmt.Printf("can't create meta dataset for %v err: %v", name, err)
	}
	defer dset.Close()

	err = dset.Write(&arr)
	if err != nil {
		fmt.Printf("can't write meta dataset for %v, err: %v", name, err)
	}
}

func saveTime(name string, grp *hdf5.Group, arr []int64) {
	nmeta := uint(len(arr))

	fildims := []uint{nmeta}
	filspace, err := hdf5.CreateSimpleDataspace(fildims, fildims)
	if err != nil {
		fmt.Println("can't create meta file dataspace", err)
	}
	defer filspace.Close()
	dset, err := grp.CreateDataset(name, hdf5.T_NATIVE_INT64, filspace)
	if err != nil {
		fmt.Printf("can't create time dataset for %v err: %v", name, err)
	}
	defer dset.Close()

	err = dset.Write(&arr)
	if err != nil {
		fmt.Printf("can't write time dataset for %v, err: %v", name, err)
	}
}

func getNodename(euiStr string, dbname string) (nodename string) {
	db, err := sql.Open("sqlite3", dbname)
	defer db.Close()
	if err != nil {
		fmt.Println("can't open db file", dbname)
		return "UNKNOWN"
	}
	rows, err := db.Query("SELECT nodename FROM nodes WHERE eui = ?;", euiStr)
	defer rows.Close()
	for rows.Next() {
		var n string
		if err := rows.Scan(&n); err != nil {
			fmt.Printf("err scanning db for %v", euiStr)
		}
		return n
	}
	return "UNKNOWN"
}

func saveAttr(name string, grp *hdf5.Group, val interface{}) (err error) {
	fmt.Printf("val %v name %v type %T\n", val, name, val)
	aspace, err := hdf5.CreateDataspace(hdf5.S_SCALAR)
	if err != nil {
		fmt.Println("can't create attr dataspace")
		return err
	}
	defer aspace.Close()

	var attr *hdf5.Attribute
	switch val.(type) {
	case float64:
		attr, err = grp.CreateAttribute(name, hdf5.T_IEEE_F64LE, aspace)
	case string:
		attr, err = grp.CreateAttribute(name, hdf5.T_GO_STRING, aspace)
	case uint:
		attr, err = grp.CreateAttribute(name, hdf5.T_NATIVE_UINT32, aspace)
	default:
		fmt.Printf("unknown datatype for attr val %v %T", val, val)
		return err
	}
	if err != nil {
		fmt.Println("can't create attr")
	}
	defer attr.Close()

	switch v := val.(type) {
	case float64:
		err = attr.Write(&v, hdf5.T_IEEE_F64LE)
	case string:
		err = attr.Write(&v, hdf5.T_GO_STRING)
	case uint:
		err = attr.Write(&v, hdf5.T_NATIVE_UINT32)
	}
	if err != nil {
		fmt.Println("can't write attr")
	}

	return nil
}
