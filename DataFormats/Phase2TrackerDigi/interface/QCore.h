#ifndef QCORE_H
#define QCORE_H
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
	int ccol_in,
	int qcrow_in,
	bool isneighbour_in,
	bool islast_in,
	std::vector<int> adcs_in
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
  std::vector<bool> toRocCoordinates(std::vector<bool>& hitmap);

  std::vector<bool> intToBinary(int num, int length);

  bool containsHit(std::vector<bool>& hitmap);

  std::vector<bool> getHitmapCode(std::vector<bool> hitmap);

};
#endif // QCORE_H
