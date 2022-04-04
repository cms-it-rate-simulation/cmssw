#include "../interface/QCore.h"
#include <cmath>
#include <vector>
#include <iostream>

//A 4x4 region of hits in sensor coordinates
/*QCore::QCore(
             int ccol_in,
             int qcrow_in,
             bool isneighbour_in,
             bool islast_in,
             std::vector<int> adcs_in
             ) :
    adcs(adcs_in),
    isneighbour(isneighbour_in),
    islast(islast_in),
    qcrow(qcrow_in),
    ccol(ccol_in)
    {
    };*/

QCore::QCore(int rocid, int ccol, int qcrow, bool isneighbour, bool islast, std::vector<int> adcs) {
  rocid_ = rocid;
  ccol_ = ccol;
  qcrow_ = qcrow;
  isneighbour_ = isneighbour;
  islast_ = islast;
  adcs_ = adcs;
}

//Returns the hitmap for the QCore in the 4x4 sensor coordinates
std::vector<bool> QCore::getHitmap() {
    	//assert(adcs.size()==16);
    
	std::vector<bool> hitmap = {};

    	//hitmap.reserve(16);
    	for(auto adc : adcs_) {
        	hitmap.push_back(adc > 0);
    	}
	//assert(hitmap.size()==16);

    	return toRocCoordinates(hitmap);
}

//Returns the bit code associated with the QCore
std::vector<bool> QCore::encodeQCore(bool is_new_col) {
	std::vector<bool> code = {};

	if(is_new_col) {
		std::vector<bool> col_code = intToBinary(ccol_, 6);
		code.insert(code.end(), col_code.begin(), col_code.end());
	}

	code.push_back(islast_);
	code.push_back(isneighbour_);

	if(!isneighbour_) {
		std::vector<bool> row_code = intToBinary(qcrow_, 8);
		code.insert(code.end(), row_code.begin(), row_code.end());
	}

	std::vector<bool> hitmap_code = getHitmapCode(getHitmap());
	code.insert(code.end(), hitmap_code.begin(), hitmap_code.end());

	for(auto adc : adcs_) {
		if(adc != 0) {
			std::vector<bool> adc_code = intToBinary(adc, 4);
			code.insert(code.end(), adc_code.begin(), adc_code.end());
		}
	}

	return code;	
}

//Takes in a hitmap in sensor corrdinates with 4x4 regions and converts it to readout chip coordinates with 2x8 regions
std::vector<bool> QCore::toRocCoordinates(std::vector<bool>& hitmap) {
	std::vector<bool> ROC_hitmap(16, 0);
	//ROC_hitmap.reserve(hitmap.size());

	for(size_t i = 0; i < hitmap.size(); i++) {
		int row = std::floor(i / 4);
		int col = i % 4;

		int new_row;
		int new_col;
		if(row % 2 == 0) {
                        new_row = row / 2;
			new_col = 2 * col;
                } else {
                        new_row = std::floor(row / 2);
			new_col = 2 * col + 1;
		}

		int new_index = 8 * new_row + new_col;
		ROC_hitmap[new_index] = hitmap[i];
	}

	return ROC_hitmap;
}

//Converts an integer into binary, and formats it with the given length
std::vector<bool> QCore::intToBinary(int num, int length) {
	int n = num;
	std::vector<bool> bi_num = {};
	//bi_num.reserve(length);
	
	for(int i = 0; i < length; i++) {
		bi_num.push_back(0);
	}

	for(int i = length; i > 0; i--) {
		if(n >= pow(2, i - 1)) {
			bi_num[length - i] = 1;
			n -= pow(2,i - 1);
		} else {
			bi_num[length - i] = 0;
		}
	}

	return bi_num;
}

//Takes in a hitmap and returns true if it contains any hits
bool QCore::containsHit(std::vector<bool>& hitmap) {
        for(size_t i = 0; i < hitmap.size(); i++) {
                if(hitmap[i]) {
                        return true;
                }
        } return false;
}

//Returns the code associated with a given hitmap
std::vector<bool> QCore::getHitmapCode(std::vector<bool> hitmap) {
        std::vector<bool> code = {};
        
	if(hitmap.size() == 1) {
                return code;
        }

        std::vector<bool> left_hitmap = std::vector<bool>(hitmap.begin(), hitmap.begin() + hitmap.size() / 2);
        std::vector<bool> right_hitmap = std::vector<bool>(hitmap.begin() + hitmap.size() / 2, hitmap.end());

        bool hit_left = containsHit(left_hitmap);
        bool hit_right = containsHit(right_hitmap);

        if(hit_left && hit_right) {
		code.push_back(1);
		code.push_back(1);

		std::vector<bool> left_code = getHitmapCode(left_hitmap);
		std::vector<bool> right_code = getHitmapCode(right_hitmap);
		code.insert(code.end(),left_code.begin(),left_code.end());
		code.insert(code.end(),right_code.begin(),right_code.end());
        } else if(hit_right) {
		code.push_back(0);
		
		std::vector<bool> right_code = getHitmapCode(right_hitmap);
		code.insert(code.end(),right_code.begin(),right_code.end());
        } else if(hit_left) {
		code.push_back(1);
		code.push_back(0);

		std::vector<bool> left_code = getHitmapCode(left_hitmap);
		code.insert(code.end(),left_code.begin(),left_code.end());
        } return code;
}
