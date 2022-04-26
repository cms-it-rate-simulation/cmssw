#ifndef DataFormats_Phase2TrackerDigi_QCore_h
#define DataFormats_Phase2TrackerDigi_QCore_h
#include<vector>

class QCore{

 private:
  std::vector<int> adcs;
  bool islast_;
  bool isneighbour_;
  int rocid_;
  int ccol_;
  int qcrow_;

 public:
  QCore(
	int rocid,
	int ccol,
	int qcrow,
	bool isneighbour,
	bool islast,
	std::vector<int> adcs
	);

  QCore() {
    rocid_ = -1;; 
    islast_ = false;
    isneighbour_ = false;
    ccol_ = -1;
    qcrow_ = -1;
  }


  void setIsLast(bool islast) {
    islast_ = islast;
  }

  bool islast() const {
    return islast_;
  }

  void setIsNeighbour(bool isneighbour) {
    isneighbour_ = isneighbour;
  }

  int rocid() const {
    return rocid_;
  }
  
  //Returns the column number of the QCore
  int ccol() const {
    return ccol_;
  }
  
  //Returns the row number of the QCore
  int qcrow() const {
    return qcrow_;
  }
  
  std::vector<bool> getHitmap();
    
  std::vector<bool> encodeQCore(bool is_new_col);
  
  const bool operator<(const QCore& other) {
    if (ccol_==other.ccol_) {
      return ccol_ < other.ccol_;
    } else {
      return qcrow_ < other.qcrow_;
    }
  }


 private:

  std::vector<int> adcs_;
  bool islast_;
  bool isneighbour_;
  int rocid_;
  int ccol_;
  int qcrow_;

  std::vector<bool> toRocCoordinates(std::vector<bool>& hitmap);

  std::vector<bool> intToBinary(int num, int length);

  bool containsHit(std::vector<bool>& hitmap);

  std::vector<bool> getHitmapCode(std::vector<bool> hitmap);

};
#endif // QCORE_H
