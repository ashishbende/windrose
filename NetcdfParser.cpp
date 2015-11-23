//============================================================================
// Name        : CMPE275Project1.cpp
// Author      : Ashish,Ashok,Nishant,Prashanth
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cstdio>
#include "math.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <cstdlib>
#include "zlib.h"
#include "zconf.h"
#include "omp.h"
#include <fstream>
#include <netcdf>
#include <list>
#include "WindroseModel.h"
#include <map>
#include <sstream>
#include <stdlib.h>
#include <dirent.h>
#include <limits>
#include <typeinfo>
#include <unordered_map>
#include <string>
#include <exception>
#include <Python.h>

using namespace std;

//using namespace netCDF;
using namespace netCDF::exceptions;
// Return this in event of a problem.
//static const int NC_ERR = 2;

//int numData = 0;
#define pi 3.14159265358979323846
#define earthRadiusKm 6371.0
#define NUM_THREADS 4
const int sectors = 16;
const int numSpds = 6;
static unsigned int NLAT = 120000;
static list<WindroseModel*> lstAttributes;
static map<string, WindroseModel> model;
static int arrayToReturn[sectors][numSpds];
static WindroseModel objToReturn;

int calcSpeedBin(float windSpd) {
	if (windSpd == std::numeric_limits<float>::max()
			|| windSpd == std::numeric_limits<float>::min()) {
		return -1;
	}

	if (windSpd >= 0 and windSpd <= 5)
		return 0;
	else if (windSpd > 5 and windSpd <= 10)
		return 1;
	else if (windSpd > 10 and windSpd <= 15)
		return 2;
	else if (windSpd > 15 and windSpd <= 20)
		return 3;
	else if (windSpd > 20 and windSpd <= 25)
		return 4;
	else
		return 5;
}

string Convert(float number) {
	ostringstream buff;
	buff << number;
	return buff.str();
}

int calcDirectionBin(float windDir) {
	// 0-360 degree - cut into a linear line 0-359
	//cout<<"Wind Direction"<<windDir<<endl;
	if (windDir == std::numeric_limits<float>::max()
			|| windDir == std::numeric_limits<float>::min())
		windDir = -1;
	float sec = -1;
	if (windDir > 360)
		return (int) sec;
	sec = (windDir / 360) * 16;
	//cout<<sec<< " " <<(int) sec;
	return (int) sec;
}

float covertToFloat(const std::string& s) {
	float f;
	std::istringstream iss(s);
	iss >> f;
	return f;
}

void printWindroseMatrix(string stationId, map<string, WindroseModel> mo) {
	WindroseModel roseModel = mo.find(stationId)->second;
	vector<float> windDirection;
	vector<float> windSpeed;
	/*istringstream f1(roseModel.WindDirection);
	string temp1;
	while (getline(f1, temp1, ',')) {
		windDirection.push_back(covertToFloat(temp1));
	}

	istringstream f2(roseModel.WindSpeed);
	string temp2;
	while (getline(f2, temp2, ',')) {
		windSpeed.push_back(covertToFloat(temp2));
	}*/

	int size = roseModel.windSpeed.size();
	//passed into fn() or ....
	//construct a wind rose data
	int wr[sectors][numSpds] = { };

	for (int i = 0; i < size; i++) {
		int s = calcSpeedBin(roseModel.windSpeed[i]);
		int d = calcDirectionBin(roseModel.windDirection[i]);
		if (s != -1 || d != -1) {
			wr[d][s]++;
		}
	}

	for (int i = 0; i < sectors; i++) {
		for (int j = 0; j < numSpds; j++) {
			cout << wr[i][j] << "\t\t";
		}
		cout << endl << endl;
	}

}

int** printWindroseMatrixHash(WindroseModel mo) {
	cout<<"Calculating windriestart"<<endl;
	WindroseModel roseModel = mo;
	vector<float> windDirection;
	vector<float> windSpeed;
	//istringstream f1(roseModel.WindDirection);
	//string temp1;
	//cout<<"Calculating windrie1"<<endl;
	/*while (getline(f1, temp1, ',')) {
		windDirection.push_back(covertToFloat(temp1));
	}*/
	//cout<<"Calculating windrie2"<<endl;

	//istringstream f2(roseModel.WindSpeed);
	//string temp2;
	/*while (getline(f2, temp2, ',')) {
		windSpeed.push_back((covertToFloat(temp2)));
	}*/

	int size = mo.windSpeed.size();
	//passed into fn() or ....
	//construct a wind rose data
	int** wr;

	wr = new int*[16];
	//cout<<"Calculating windrie3"<<endl;
	for (int i = 0; i < 16; i++) {
		wr[i] = new int[6];
		for (int j = 0; j < 6; j++) {
			wr[i][j] = 0;
		}
	}
	//cout<<"Calculating windrie4"<<endl;
	for (int i = 0; i < size; i++) {
		//cout<<i<<","<<size<<endl;
		int s = calcSpeedBin(mo.windSpeed[i]);

		int d = calcDirectionBin(mo.windDirection[i]);
		if (s != -1 && d != -1 && mo.windSpeed[i] != 0.0f) {
			wr[d][s] = wr[d][s] + 1;

		}
	}
	//cout<<"Calculating windrie5"<<endl;
	return wr;

}

bool gzipInflate(const std::string& compressedBytes,
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
bool loadBinaryFile(const std::string& filename, std::string& contents) {
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

vector<string> getAllFilesinDirectory(char* directoryPath) {
	DIR *pDIR;
	vector<string> files;
	struct dirent *entry;
	if (pDIR = opendir(directoryPath)) {
		while (entry = readdir(pDIR)) {
			if (strcmp(entry->d_name, ".") != 0
					&& strcmp(entry->d_name, "..") != 0)
				files.push_back(entry->d_name);
		}

		closedir(pDIR);
	}

	return files;
}

double deg2rad(double deg) {
	return (deg * pi / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
	return (rad * 180 / pi);
}

double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
	//printf("Calculating0");
	double lat1r, lon1r, lat2r, lon2r, u, v;
	//cout << "Calculating1";
	lat1r = deg2rad(lat1d);
	//cout << "Calculating1";
	lon1r = deg2rad(lon1d);
	//cout << "Calculating1";
	lat2r = deg2rad(lat2d);
	//cout << "Calculating1";
	lon2r = deg2rad(lon2d);
	//cout << "Calculating1";
	u = sin((lat2r - lat1r) / 2);
	//cout << "Calculating1";
	v = sin((lon2r - lon1r) / 2);
	//cout << "Calculating1";
	return 2.0 * earthRadiusKm
			* asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

string calculateAverageSpeed(int sector, int** windRoseMatrix) {
	float intAvgSpeed;
	float totalRecords;
	float lowerLimit = 0.0;
	float higherLimit = 0.0;
	float intOccurences = 0.0;

	for (int i = 0; i < 5; i++) {
		switch (i) {
		case 0:
			lowerLimit += (windRoseMatrix[sector][i] * 0);
			higherLimit += (windRoseMatrix[sector][i] * 5);
			intOccurences += windRoseMatrix[sector][i];
			break;
		case 1:
			lowerLimit += (windRoseMatrix[sector][i] * 6);
			higherLimit += (windRoseMatrix[sector][i] * 10);
			intOccurences += windRoseMatrix[sector][i];
			break;
		case 2:
			lowerLimit += (windRoseMatrix[sector][i] * 11);
			higherLimit += (windRoseMatrix[sector][i] * 15);
			intOccurences += windRoseMatrix[sector][i];
			break;
		case 3:
			lowerLimit += (windRoseMatrix[sector][i] * 16);
			higherLimit += (windRoseMatrix[sector][i] * 25);
			intOccurences += windRoseMatrix[sector][i];
			break;
		case 4:
			lowerLimit += (windRoseMatrix[sector][i] * 26);
			higherLimit += (windRoseMatrix[sector][i] * 26);
			intOccurences += windRoseMatrix[sector][i];
			break;
		}
		//lowerLimit = windRoseMatrix[sector][i] ;

		if (intOccurences == 0) {
			return "0,0";
		} else {
			return Convert(lowerLimit / intOccurences) + ","
					+ Convert(higherLimit / intOccurences);
		}
	}
}

std::string calulcateSpeedFromOriginToDestination(WindroseModel o,
		WindroseModel d) {
	cout << "Calculating windrie0.0" << endl;
	int** windRoseMatrix = printWindroseMatrixHash(o);
	cout << "windrie completed" << endl;
	//arrayToReturn = windRoseMatrix;
	float longDist = o.Longitude - d.Longitude;
	float latDist = o.Latitude - d.Latitude;
	float angle = atan2(longDist, latDist) * 180 / pi;
	cout << "Angle between 2 points is = " << angle;
	if (angle < 0.0) {
		angle += 360.0;
	}
	angle = calcDirectionBin(angle);
	cout << "Angle secotor " << angle << endl;
	angle = 6;
	string avg = calculateAverageSpeed(angle, windRoseMatrix);
	vector<float> avgSpeed;
	istringstream f1(avg);
	string temp1;
	//cout<<"Calculating windrie1"<<endl;
	while (getline(f1, temp1, ',')) {
		avgSpeed.push_back(covertToFloat(temp1));
	}
	cout << "Lower Limit " << avgSpeed[0] << endl;
	cout << "higher Limit " << avgSpeed[1] << endl;
	//cout<<"higher Limit "<<avgSpeed[1]<<endl;
	cout << "Calculating Average :" << endl;
	float average1;
	float average2;
	//average1 = avgSpeed[0] + avgSpeed[1];
	//cout<<"Addition result " <<average <<endl;
	//average = average/2.0;
	//cout<<"Division result " <<average <<endl;

	cout << "Calculating distance of earth" << endl;
	cout << "Calculating distance of earth" << endl;
	cout << "Calculating distance of earth" << endl;

	double distance = 0;
	distance = distanceEarth((double) o.Latitude, (double) o.Longitude,
			(double) d.Latitude, (double) d.Longitude);
	cout << "Distance of earth " << distance << endl;
	/*average1 = (float)(avgSpeed[0]);
	 cout<<"Average speed o = "<<average1;
	 if (avgSpeed[0] > 0) {
	 cout<<"In calculating Avg speed 1 ";
	 average1 = avgSpeed[0] * 3.6;
	 double time = distance / average1;
	 cout << "It will reach in these minimum hours " << time;
	 } else {
	 cout << "It Never reach in miminum time";
	 }
	 if (avgSpeed[1] > 0.0) {
	 cout<<"In calculating Avg speed 2 ";
	 average2 = avgSpeed[1] * 3.6;
	 double time = distance / average2;
	 cout << "It will reach in these maximum hours " << time;
	 } else {
	 cout << "It Never reach in maximum time";
	 }*/
	cout << "Calculating distance of earth" << endl;
	return "";
}
void printTime(std::string strMessage)
{
	time_t rawtime;
	struct tm * timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	cout<<strMessage << ":"<<asctime(timeinfo)<<endl;

}
WindroseModel mainProgram(int argc, char **argv) {
	time_t rawtime;
	struct tm * timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current start time: %s", asctime(timeinfo));

	//char* directoryPath = argv[1];
	char* directoryPath = "/home/nishant/windrosedata/data";
	float lat = 60.4522;
	float lng = 5.135;
	std::unordered_map<std::string, WindroseModel> um[NUM_THREADS];
	vector<string> files = getAllFilesinDirectory(directoryPath);
	int fileSize = files.size();
	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel
	{

		int id = omp_get_thread_num();
		cout << "Starting thread - " << id;
		for (int i = id; i < fileSize; i += NUM_THREADS) {
			char* directoryPath2 = "/home/nishant/windrosedata/data";
			std::string fileData;
			std::string fileName(directoryPath2);
			fileName = fileName + "/" + files[i];
			int randomNumber = rand() % 999999;
			stringstream ss1;
			ss1 << randomNumber;
			std::string outputFilePath = "/dev/shm/" + ss1.str();
			cout << outputFilePath << endl;


			if (!loadBinaryFile(fileName, fileData)) {
				printf("Error loading input file.");
				//return 0;
			}
			/*time ( &rawtime );
			 timeinfo = localtime ( &rawtime );
			 printf ( "Starting inflating  file : %s %d" , asctime (timeinfo),id);
			 cout<<endl;*/
			std::string data;
			if (!gzipInflate(fileData, data)) {
				printf("Error decompressing file.");
				//return 0;
			}
			/*time ( &rawtime );
			 timeinfo = localtime ( &rawtime );
			 printf ( "inflation complete : %s %d" , asctime (timeinfo),id);
			 cout<<endl;*/
			std::ofstream out(outputFilePath.c_str());
			out << data;
			out.close();

			cout << "Opening file in thread " << id<< endl;
			//#pragma omp barrier

			//#pragma omp barrier

			float index[NLAT];
			float index2[NLAT];
			char index3[NLAT][6];
			float index6[NLAT];
			float index7[NLAT];
			char index4[NLAT][11];
			char index5[NLAT][12];

//#pragma omp critical
			//{
				netCDF::NcVar windDir, windSpeed, stationId, stationType,
						providerId, latitude, longitude;
				//cout << "Read Char for thread- "<<id << endl;

				/*time ( &rawtime );
				 timeinfo = localtime ( &rawtime );
				 printf ( "Reading in thread start : %s %d" , asctime (timeinfo),id);*/
				netCDF::NcFile dataFile = netCDF::NcFile(outputFilePath,
						netCDF::NcFile::read);
				windDir = dataFile.getVar("windDir");
				windSpeed = dataFile.getVar("windSpeed");
				stationId = dataFile.getVar("stationId");
				latitude = dataFile.getVar("latitude");
				longitude = dataFile.getVar("longitude");
				//stationType = dataFile.getVar("stationType");
				//providerId = dataFile.getVar("providerId");

				//windSpeed.getVar()

				windSpeed.getVar(index);

				windDir.getVar(index2);

				stationId.getVar(index3);
				//cout << "Read Char2.3" << endl;

				//stationType.getVar(index4);
				//cout << "Read Char2" << endl;
				//providerId.getVar(index5);
				latitude.getVar(index6);
				longitude.getVar(index7);
				/*time ( &rawtime );
				 timeinfo = localtime ( &rawtime );
				 printf ( "Reading Complete in thread start : %s %d" , asctime (timeinfo),id);*/
			//}

			for (int i = 0; i < NLAT; i++) {
				//cout<<"reading "<< i<<endl ;
				std::string str = index3[i];
				string str1(index4[i]);
				/*if (str != "\0" && str != "") {
				 oCSV( 0, i+1 ) = str;
				 oCSV( 1, i+1 ) = Convert(index[i]);
				 oCSV( 2, i+1 ) = Convert(index2[i]);
				 oCSV( 3, i+1 ) = index4[i];
				 oCSV( 4, i+1 ) = index5[i];
				 oCSV( 5, i+1 ) = Convert(index6[i]);
				 oCSV( 6, i+1 ) = Convert(index7[i]);
				 }*/

					if (um[id].find(str) == um[id].end()) {
						//um.insert
						if(index[i]>0.0)
						{
							//string windS = std::to_string(calcSpeedBin(index[i]));
						um[id][str] = WindroseModel(index6[i], index7[i]);
						um[id][str].windSpeed.push_back(index[i] * 3.6);
						um[id][str].windDirection.push_back(index2[i]);
						}
						//cout<<"New object created"
						//um.insert(str,WindroseModel(Convert(index[i]),Convert(index2[i]),index6[i],index7[i]));
					} else {
						//cout<<"New object Not created"<<endl;
						if(index[i]>0.0)
						{
							//string windS = std::to_string(calcSpeedBin(index[i]));
							WindroseModel windRMD = um[id][str];
							/*windRMD.WindDirection = windRMD.WindDirection + ","
									+  Convert(index2[i]);
							windRMD.WindSpeed = windRMD.WindSpeed + ","
									+ std::to_string((index[i]* 3.6));*/
							um[id][str] = windRMD;


							um[id][str].windSpeed.push_back(index[i] * 3.6);
							um[id][str].windDirection.push_back(index2[i]);
						}

				}

			}

			//oCSV.save(csvFile.c_str());
			/*map<string, WindroseModel> mod2; //= getStationsFromCSV("AR173", csvFile);
			 map<string, WindroseModel> mod; //= getStationsFromCSV("AS070", csvFile);
			 mod2["AR173"] = um[0]["AR173"];
			 mod["AS070"] = um[0]["AS070"];
			 cout << endl << "Output WindSpeed of AS070 " << mod["AS070"].WindSpeed
			 << endl;
			 cout << endl << "Output WindDirec of AS070 "
			 << mod["AS070"].WindDirection << endl;

			 cout << endl << "Output WindSpeed of AR173 " << mod2["AR173"].WindSpeed
			 << endl;
			 cout << endl << "Output WindDirec of AR173 "
			 << mod2["AR173"].WindDirection << endl;

			 //printWindroseMatrix("AS070", mod);
			 cout << mod2["AR173"].Latitude << "," << mod2["AR173"].Longitude
			 << endl;
			 *//*cout<< mod["AS070"].Latitude<<","<<mod["AR120"].Longitude<<endl; //39.9253,-81.435
			 cout<< mod["AS070"].WindDirection<<endl;
			 cout<< mod["AS070"].WindSpeed<<endl;*/
			//cout<<endl<<endl;
			//cout<<"SIZE OF MAP IS "<<model.size()<<endl;
			//cout<<"MAX SIZE OF MAP "<<model.max_size()<<endl;
			//float longDist = mod2["AR173"].Longitude - mod["AS070"].Longitude;
			//float latDist= mod2["AR173"].Latitude - mod["AS070"].Latitude;
			if (remove(outputFilePath.c_str()) != 0)
				perror("Error deleting file");
			else
				puts("File successfully deleted");
		}
		//#pragma omp barrier

	}
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current End time: %s", asctime(timeinfo));
	int option = 1;
	while (option == 1) {
		std: string station1, station2;
		cout << "Enter Station id for origin" << endl;
		//std::cin>>station1;
		station1 = "AR173";
		cout << "Enter Station id for destination" << endl;
		//std::cin>>station2;
		station2 = "AS070";
		WindroseModel o = WindroseModel();
		WindroseModel d = WindroseModel();
		;
		for (int i = 0; i < NUM_THREADS; i++) {
			//cout<<"Wind speed of station in thread " << i<<" is " <<um[i][station2].WindSpeed<<endl;
			if (um[i].find(station1) != um[i].end()) {
				if (o.windDirection.size() == 0) {
					o = WindroseModel(um[i][station1].Latitude,
							um[i][station1].Longitude);
					o.windSpeed.swap(um[i][station1].windSpeed);
					o.windDirection.swap(um[i][station1].windDirection);
				} else {
					/*o.WindDirection = um[i][station1].WindDirection + ","
							+ o.WindDirection;
					o.WindSpeed = um[i][station1].WindSpeed + "," + o.WindSpeed;*/
					o.windSpeed.insert(o.windSpeed.end(),um[i][station1].windSpeed.begin(),um[i][station1].windSpeed.end());
					o.windDirection.insert(o.windDirection.end(),um[i][station1].windDirection.begin(),um[i][station1].windDirection.end());
				}
			}
			if (um[i].find(station2) != um[i].end()) {
				if (d.windDirection.size() == 0) {
					d = WindroseModel(um[i][station2].Latitude,
							um[i][station2].Longitude);
					d.windSpeed.swap(um[i][station2].windSpeed);
					d.windDirection.swap(um[i][station2].windDirection);
				}

				else {
					/*d.WindDirection = um[i][station2].WindDirection + ","
							+ d.WindDirection;
					d.WindSpeed = um[i][station2].WindSpeed + "," + d.WindSpeed;*/
					d.windSpeed.insert(d.windSpeed.end(),um[i][station2].windSpeed.begin(),um[i][station2].windSpeed.end());
					d.windDirection.insert(d.windDirection.end(),um[i][station2].windDirection.begin(),um[i][station2].windDirection.end());
				}
			}
		}
		objToReturn = d;
		/*cout << "Wind Direction for origin " << o.WindDirection << endl;
		cout << "Wind Speed for origin " << o.WindSpeed << endl;
		cout << "Wind Direction for destination " << d.WindDirection << endl;
		cout << "Wind Speed for destination " << d.WindSpeed << endl;*/
		string str = calulcateSpeedFromOriginToDestination(d, o);
		//cout<<"Enter 0 to exit";
		//cin>>option;
		return objToReturn;
	}
}

int main(int argc,char** argv)
{
	WindroseModel obj = mainProgram(argc, argv);

	PyObject* py_imp_str;
	PyObject* py_imp_handle;
	PyObject* py_imp_dict; //borrowed
	PyObject* py_imp_load_source; //borrowed
	PyObject* py_dir; //stolen
	PyObject* py_lib_name; //stolen
	PyObject* py_args_tuple;
	PyObject* py_lib_mod;
	PyObject* py_func;
	PyObject* py_ret;
	PyObject* py_lib_mod_dict;
	PyObject *expr1,*expr2, *args1;

	    wstring filestr(L"/home/nishant/workspacecpp/netcdf/src/windrose.py");
	    wstring functionstr(L"plotit2");

	    const wchar_t * file = filestr.c_str();
	    const wchar_t * function = functionstr.c_str();

	    Py_Initialize();

	    py_dir = PyUnicode_FromWideChar(file, wcslen(file) );
	    py_imp_str = PyString_FromString("imp");
	    py_imp_handle = PyImport_Import(py_imp_str); //normal python import for imp
	    if(!py_imp_handle)
	    {
	    	cout<<"Error";
	    }
	    py_imp_dict = PyModule_GetDict(py_imp_handle); //borrowed
	    py_imp_load_source = PyDict_GetItemString(py_imp_dict, "load_source"); //borrowed
	    py_lib_name = PyUnicode_FromWideChar(function, wcslen(function));

	    py_args_tuple = PyTuple_New(2);
	    PyTuple_SetItem(py_args_tuple, 0, py_lib_name); //stolen
	    PyTuple_SetItem(py_args_tuple, 1, py_dir); //stolen
	    args1 = PyTuple_New(2);
	    py_lib_mod = PyObject_CallObject(py_imp_load_source, py_args_tuple);
	    py_lib_mod_dict = PyModule_GetDict(py_lib_mod); //borrowed



	    printTime("Starting Iteration at time ");
		ostringstream oss;
		copy(obj.windSpeed.begin(), obj.windSpeed.end()-1,ostream_iterator<float>(oss, ","));
		oss << obj.windSpeed.back();

		std::ostringstream oss2;
		copy(obj.windDirection.begin(), obj.windDirection.end()-1,ostream_iterator<float>(oss2, ","));
		oss2 << obj.windDirection.back();

		std::string str1 = oss.str() ;
		std::string str2 = oss2.str() ;

		cout<<str1<<endl;
		cout<<str2<<endl;
		printTime("Ending Iteration at time ");


	    char *a1 = new char[str1.size() + 1];
		a1[str1.size()] = 0;
		memcpy(a1, str1.c_str(), str1.size());
		char *a2 = new char[str2.size() + 1];
		a2[str2.size()] = 0;
		memcpy(a2, str2.c_str(), str2.size());

		expr1 = PyString_FromString(a1);
		expr2 = PyString_FromString(a2);

		PyTuple_SetItem(args1, 0, expr1);
		PyTuple_SetItem(args1, 1, expr2);



	    py_func = PyDict_GetItem(py_lib_mod_dict, py_lib_name);
	    py_ret = PyObject_CallFunction(py_func,"ss",a1,a2);




	    if (py_ret)
	    {
	        printf("PyObject_CallFunction from wmain was successful!\n");
	        Py_DECREF(py_ret);
	    }
	    else
	        printf("PyObject_CallFunction from wmain failed!\n");

	    Py_DECREF(py_imp_str);
	    Py_DECREF(py_imp_handle);
	    Py_DECREF(py_args_tuple);
	    Py_DECREF(py_lib_mod);
	    Py_Finalize();

	    fflush(stderr);
	    fflush(stdout);
	    return 0;
}
