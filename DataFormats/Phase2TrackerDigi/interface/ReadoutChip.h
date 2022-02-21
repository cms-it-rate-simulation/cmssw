#ifndef READOUTCHIP_H
#define READOUTCHIP_H

#include <vector>
#include <utility>
#include <string>
#include "QCore.h"
#include "Hit.h"

class ReadoutChip {
	std::vector<Hit> hitList;
	int rocnum_;

public:
	ReadoutChip(int rocnum, std::vector<Hit> hl);

	//Returns the number of hits on the readout chip
	unsigned int size();

	int rocnum() const {
	  return rocnum_;
	}

	std::vector<QCore> getOrganizedQCores();

	std::vector<bool> getChipCode();

private:
	std::pair<int,int> get_QCore_pos(Hit hit);

	QCore getQCoreFromHit(Hit pixel);

	std::vector<QCore> rem_duplicates(std::vector<QCore> qcores);

	std::vector<QCore> organize_QCores(std::vector<QCore> qcores);
};

#endif
