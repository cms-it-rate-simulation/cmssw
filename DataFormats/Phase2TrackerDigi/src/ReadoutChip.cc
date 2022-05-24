#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include "../interface/QCore.h"
#include "../interface/ReadoutChip.h"
#include "../interface/Hit.h"

bool ReadoutChip::endOfStreamMarker = true;

ReadoutChip::ReadoutChip(int rocnum, std::vector<Hit> hitList) {
  hitList_ = hitList;
  rocnum_ = rocnum; 
}

unsigned int ReadoutChip::size() {
	return hitList_.size();
}

//Takes in list of hits and organizes them into the 4x4 QCores that contains them
std::vector<QCore> ReadoutChip::getOrganizedQCores() {
  std::cout << "In getOrganizedQCores" <<std::endl;
        std::vector<QCore> qcores = {};
	bool qcore_already_exists;
	std::pair<int,int> qcore_pos;

        for(const auto& hit : hitList_) {
		qcore_already_exists = false;
		qcore_pos = getQCorePos(hit);

		for(size_t i = 0; i < qcores.size(); i++) {
			if(qcores[i].qcrow() == qcore_pos.first && qcores[i].ccol() == qcore_pos.second) {
				qcore_already_exists = true;
			}
		}

		if(!qcore_already_exists) {
                	qcores.push_back(getQCoreFromHit(hit));
		}
        }

        return linkQCores(organizeQCores(qcores));
}

//Returns the encoding of the readout chip
std::vector<bool> ReadoutChip::getChipCode(int event, bool aurora) {
        std::vector<bool> code = intToBinary(event, 8);

        if(hitList_.size() > 0) {
                std::vector<QCore> qcores = getOrganizedQCores();

		bool is_new_col = true;

                for(auto& qcore : qcores) {
                	std::vector<bool> qcore_code = qcore.encodeQCore(is_new_col);
			code.insert(code.end(), qcore_code.begin(), qcore_code.end());
			
			is_new_col = qcore.islast();
                }
        }

	if(aurora) {
		auroraFormat(code);
	}

        return code;
}

//Returns the position (row,col) of the 4x4 QCore that contains a given hit
std::pair<int,int> ReadoutChip::getQCorePos(Hit hit) {

        int row = hit.row() / 4;
        int col = hit.col() / 4;

        return {row,col};
}

//Takes in a hit and returns the 4x4 QCore that hit is a part of
QCore ReadoutChip::getQCoreFromHit(Hit pixel) {
        std::vector<int> adcs =
                {0,0,0,0,
		0,0,0,0,
                0,0,0,0,
		0,0,0,0};

        std::pair<int,int> pos = getQCorePos(pixel);

        for(const auto& hit : hitList_) {
                if(getQCorePos(hit) == pos) {
                        int i = (4 * (hit.row() % 4) + (hit.col() % 4) + 8) % 16;
                        adcs[i] = hit.adc();
                }
        }

        QCore qcore(0, pos.second, pos.first, false, false, adcs);

        return qcore;
}

//Returns a list of the qcores with hits arranged by increasing column then row numbers
std::vector<QCore> ReadoutChip::organizeQCores(std::vector<QCore> qcores) {
        std::vector<QCore> organized_list = {};

        while(qcores.size() > 0) {
                int min = 0;

                for(size_t i = 1; i < qcores.size(); i++) {
                        if(qcores[i].ccol() < qcores[min].ccol()) {
                                min = i;
                        } else if(qcores[i].ccol() == qcores[min].ccol() && qcores[i].qcrow() < qcores[min].qcrow()) {
                                min = i;
                        }
                }

                organized_list.push_back(qcores[min]);
                qcores.erase(qcores.begin() + min);
        }

        return organized_list;
}

//Takes in an oranized list of qcores and sets the islast and isneighbor properties of those qcores
std::vector<QCore> ReadoutChip::linkQCores(std::vector<QCore> qcores) {
  std::cout << "In link_QCores size " << qcores.size() << std::endl;
	for(size_t i = 1; i < qcores.size(); i++) {
		if(qcores[i].qcrow() == qcores[i - 1].qcrow()) {
			qcores[i].setIsNeighbour(true);
		}
	}

	std::cout << "Here001" << std::endl;

	if (qcores.size()>0) {

	  //size is unsigned so if size is zero size()-1 is a huge number...
	  //Hence this needs to be procted
	  for(size_t i = 0; i < qcores.size() - 1; i++) {
	    if(qcores[i].ccol() != qcores[i + 1].ccol()) {
	      qcores[i].setIsLast(true);
	    }
	  }
	  
	  std::cout << "Here002" << std::endl;
	  
	  qcores[qcores.size() - 1].setIsLast(true);
	}

	std::cout << "Here003" << std::endl;

	return qcores;
}

//Converts an integer into binary, and formats it with the given length
std::vector<bool> ReadoutChip::intToBinary(int num, int length) {
	int n = num;
	std::vector<bool> bi_num = {};
	
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

//Takes in an unformatted encoding and adds aurora formatting to it
void ReadoutChip::auroraFormat(std::vector<bool>& code) {
	addOrphanBits(code);
}

//Takes in an encoding and makes its length equal to the nearest multiple of 64
//above its current length by adding 0's to the end of the code
void ReadoutChip::addOrphanBits(std::vector<bool>& code) {
	int trailingZeros = 0;

	while(code.size() % 64 != 0) {
		code.push_back(0);
		trailingZeros++;
	}

	if(endOfStreamMarker && trailingZeros < 6) {
		for(int i = 0; i < 64; i++) {
			code.push_back(0);
		}
	}
}
