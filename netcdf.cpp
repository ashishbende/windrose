/*
 #include <ncFile.h>
 #include <ncVar.h>
 #include <ncCheck.h>
 #include <ncDim.h>
 #include <ncException.h>
 */

//#include <cstring>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include <cstdlib>
#include "zlib.h"
#include "zconf.h"
#include <fstream>
#include <netcdf>
#include "stdlib.h"
#include <map>
#include <limits>
#include <list>


using namespace std;
//using namespace netCDF;
using namespace netCDF::exceptions;
// Return this in event of a problem.
static const int NC_ERR = 2;
static long unsigned int NLAT = 300000;

class stationData
{
public:
	std::list<float> lstWindSpeed;
	std::list<float> lstWindDir;
};

bool gzipInflate1(const std::string& compressedBytes,
		std::string& uncompressedBytes) {
	if (compressedBytes.size() == 0) {
		uncompressedBytes = compressedBytes;
		return true;
	}

	uncompressedBytes.clear();

	unsigned full_length = compressedBytes.size();
	unsigned half_length = compressedBytes.size() / 2;

	unsigned uncompLength = full_length;
	char* uncomp = (char*) calloc(sizeof(char), uncompLength);

	z_stream strm;
	strm.next_in = (Bytef *) compressedBytes.c_str();
	strm.avail_in = compressedBytes.size();
	strm.total_out = 0;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;

	bool done = false;

	if (inflateInit2(&strm, (16+MAX_WBITS)) != Z_OK) {
		free(uncomp);
		return false;
	}

	while (!done) {
		// If our output buffer is too small
		if (strm.total_out >= uncompLength) {
			// Increase size of output buffer
			char* uncomp2 = (char*) calloc(sizeof(char),
					uncompLength + half_length);
			memcpy(uncomp2, uncomp, uncompLength);
			uncompLength += half_length;
			free(uncomp);
			uncomp = uncomp2;
		}

		strm.next_out = (Bytef *) (uncomp + strm.total_out);
		strm.avail_out = uncompLength - strm.total_out;

		// Inflate another chunk.
		int err = inflate(&strm, Z_SYNC_FLUSH);
		if (err == Z_STREAM_END)
			done = true;
		else if (err != Z_OK) {
			break;
		}
	}

	if (inflateEnd(&strm) != Z_OK) {
		free(uncomp);
		return false;
	}

	for (size_t i = 0; i < strm.total_out; ++i) {
		uncompressedBytes += uncomp[i];
	}
	free(uncomp);
	return true;
}

/* Reads a file into memory. */
bool loadBinaryFile1(const std::string& filename, std::string& contents) {
	// Open the gzip file in binary mode
	FILE* f = fopen(filename.c_str(), "rb");
	if (f == NULL)
		return false;

	// Clear existing bytes in output vector
	contents.clear();

	// Read all the bytes in the file
	int c = fgetc(f);
	while (c != EOF) {
		contents += (char) c;
		c = fgetc(f);
	}
	fclose(f);

	return true;
}
int getSector(float dir) {
	if(dir == std::numeric_limits<float>::max() || dir == std::numeric_limits<float>::min())
		dir = -1;
	float sec = (dir / 360) * 16;
	//cout<<sec<< " " <<(int) sec;
	return (int) sec;
}

int getSpeedRange(float spd) {
	//cout<<endl<<spd<< "  speed "<<endl;

	if(spd == std::numeric_limits<float>::max() || spd == std::numeric_limits<float>::min())
	{
		return -1;
	}
	if (spd > 0 && spd <= 5) {
		return 0;
	} else if (spd > 5 && spd <= 10) {
		return 1;
	} else if (spd > 10 && spd <= 20) {
		return 2;
	} else if (spd > 20 && spd <= 25) {
		return 3;
	} else {
		return 4;
	}

}
int main2() {
	//NcFile::FileFormat format = NcFile::FileFormat::nc4;
	//NcFile::FileMode =  NcFile::read;

	std::string fileData;
	std::string fileName;
	std::string strFileName = "";
	std::string strNewFileName = "";
	std::string strDay = "";
	std::string strSource = "/home/nishant/windrosedata/data/";
	std::string strDestination = "/home/nishant/windrosedata/extracted/";
	for (int i = 1; i < 2; i++) {
		std::ostringstream ss;
		std::string s = "";

		if (i < 10) {
			std::ostringstream ss;
			long num = i;
			ss << num;
			strDay = "0" + ss.str();
		} else {
			std::ostringstream ss;
			long num = i;
			ss << num;
			strDay = ss.str();
		}

		//strDay = ss.str();
		strFileName = "201401" + strDay;
		for (int j = 0; j < 7; j++) {
			if (j < 10) {
				std::ostringstream ss2;
				ss2 << j;
				strNewFileName = strFileName + "_0" + ss2.str() + "00";
//				/cout<<strNewFileName<<endl;
			} else {
				std::ostringstream ss2;
				ss2 << j;
				strNewFileName = strFileName + "_" + ss2.str() + "00";

			}
			strNewFileName = strSource + "01-" + strDay + "/" + strNewFileName
					+ ".gz";
			cout << strNewFileName << endl;
		}
	}
	/*if(true)
	 {
	 return 0;
	 }*/

	fileName =
			"/home/nishant/Documents/SJSU/CMPE275/Project/test/20140101_0000.gz";

	if (!loadBinaryFile1(fileName, fileData)) {
		printf("Error loading input file.");
		return 0;
	}
	std::string data;
	if (!gzipInflate1(fileData, data)) {
		printf("Error decompressing file.");
		return 0;
	}
	std::ofstream out(
			"/home/nishant/Documents/SJSU/CMPE275/Project/test/output");
	out << data;
	out.close();

	cout << "Opening file" << endl;
	//NcFile dataFile("/home/nishant/Documents/SJSU/CMPE275/Project/test/20140101_0000", NcFile::read);
	netCDF::NcFile dataFile(
			"/home/nishant/Documents/SJSU/CMPE275/Project/test/output",
			netCDF::NcFile::read);
	//cout<<"File Opened and variable count"<<dataFile.get_var("windDir")<<endl;

	netCDF::NcVar windDir, windSpeed, stationId,stationName;
	windDir = dataFile.getVar("windDir");
	windSpeed = dataFile.getVar("windSpeed");
	stationId = dataFile.getVar("stationId");
	stationName =dataFile.getVar("stationName");
	//windDir.
	//std::vector<netCDF::NcDim> dimsWindDir = windDir.getDims();
	//std::vector<netCDF::NcDim> dimsWindSpeed = windSpeed.getDims();

	float *dataValuesWindDir;
	float* dataValuesWindSpeed = new float();
	float index[NLAT];
	float index2[NLAT];
	char index3[NLAT][6];
	char index4[NLAT][51];

	//std::map<std::string, int[]> m;
	std::map<std::string, stationData> m;

	cout << "there are " << dataFile.getVarCount() << " Variables" << endl;
	cout << "there are " << dataFile.getAttCount() << " attributes" << endl;
	cout << "there are " << dataFile.getDimCount() << " dimensions" << endl;
	cout << "there are " << dataFile.getGroupCount() << " groups" << endl;
	cout << "there are " << dataFile.getTypeCount() << " types" << endl;
	//windSpeed.getVar()
	windSpeed.getVar(index);
	windDir.getVar(index2);
	stationId.getVar(index3);
	stationName.getVar(index4);
	for (int i = 0; i < NLAT; i++) {

		if (index3[i][0] != '\0') {

			//cout << i << " " << index3[i] << endl;
			std::string ind = index3[i];

			if (m.find(ind) == m.end()) {
				stationData std;
				/*for(int l=0;l<16;l++)
				{
					arr[i] = new int[5];
					for(int h=0;h<5;h++)
					{
						arr[i][h] = 4;
					}

					//cout<<arr[i][3];
				}*/
				//cout<<endl<<"Array Address -> "<<&arr;
				m.insert(std::pair<std::string, stationData>(ind, std));
			}
			//getSector(index2[i]);
			//cout<<index[i]<<" index i";
			//cout<<getSector(index2[i]) << ", " <<getSpeedRange(index[i]);
			//cout<< *(m[ind] + 5 * getSector(index2[i]) + getSpeedRange(index[i]));
			//cout<<(*m.find(ind))[getSector(index2[i])][getSpeedRange(index[i])];

			int sector = getSector(index2[i]);
			int speedrange = getSpeedRange(index[i]);
			if(sector != -1 && speedrange != -1)
			{
				//cout<<*(m[ind] + 5 * sector + speedrange);


				//int val[16][5];
				//val[0][0] = (*(m[ind]))[0][0];
				//cout<<(*(m[ind])[2][3])<<endl;
				stationData g = m[ind];
				g.lstWindDir.insert(g.lstWindDir.end(),index2[i]);
				g.lstWindSpeed.insert(g.lstWindSpeed.end(),index[i]);

				int count =0;
				for (std::list<float>::iterator it = g.lstWindDir.begin(); it != g.lstWindDir.end(); it++)
				{count++;
				    std::cout << *it << ' ';
				}
				cout<<count<<endl;

				//cout<<"Sector "<<sector<<", Speed Range "<<speedrange<<", Value "<<g[16*sector + speedrange]<<endl;

				//*(m[ind])
				//val = m[ind] ;

				//cout<< *(val+speedrange)<<endl;
				//*(m[ind] + (5 * sector) + speedrange) = ++val;
			}

			//m[ind][getSector(index2[i])][getSpeedRange(index[i])] = *m[ind][getSector(index2[i])][getSpeedRange(index[i])] + 1;

		}
	}
	cout << "there are " << dataFile.getTypeCount() << " types" << endl;
	if (remove("/home/nishant/Documents/SJSU/CMPE275/Project/test/output") != 0)
		perror("Error deleting file");
	else
		puts("File successfully deleted");
	return 0;

}

