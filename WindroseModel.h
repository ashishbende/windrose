/*
 * WindroseModel.h
 *
 *  Created on: Oct 8, 2015
 *      Author: ashok
 */

#include<string>
#include<iostream>


#ifndef WINDROSEMODEL_H_
#define WINDROSEMODEL_H_

using namespace std;

class WindroseModel{
public:
	//string WindSpeed;
	//string WindDirection;
	float Latitude;
	float Longitude;
	vector<float> windDirection;
	vector<float> windSpeed;
public:
		WindroseModel(float latitude,float longitude){
			//WindSpeed = windSpeed;
			//WindDirection = windDirection;
			Latitude = latitude;
			Longitude = longitude;
			//vector<vector<int>> a;

		}

		WindroseModel(){
			//WindSpeed = "";
			//WindDirection = "";
			Latitude = 0.0;
			Longitude = 0.0;
		}
};



#endif /* WINDROSEMODEL_H_ */
