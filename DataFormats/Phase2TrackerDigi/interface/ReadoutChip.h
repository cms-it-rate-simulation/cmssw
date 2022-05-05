#ifndef DataFormats_Phase2TrackerDigi_ReadoutChip_h
#define DataFormats_Phase2TrackerDigi_ReadoutChip_h

#include <vector>
#include <utility>
#include <string>
#include "QCore.h"
#include "Hit.h"

class ReadoutChip {

public:
	ReadoutChip(int rocnum, std::vector<Hit> hitList);

	//Returns the number of hits on the readout chip
	unsigned int size();

	int rocnum() const {
	  return rocnum_;
	}

	std::vector<QCore> getOrganizedQCores();

	std::vector<bool> getChipCode();

private:
	std::vector<Hit> hitList_;
	int rocnum_;

	std::pair<int,int> getQCorePos(Hit hit);

	QCore getQCoreFromHit(Hit pixel);

	//std::vector<bool> getChipCode();

	std::vector<QCore> remDuplicates(std::vector<QCore> qcores);

	std::vector<QCore> organizeQCores(std::vector<QCore> qcores);

	std::vector<QCore> linkQCores(std::vector<QCore> qcores);
};

#endif
