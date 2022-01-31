#ifndef QCORE_H
#define QCORE_H
#include<vector>

class QCore{

 private:
  std::vector<int> adcs;
  bool islast_;
  bool isneighbour_;
  int rocid_;
  int ccol;
  int qcrow;

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
    ccol = -1;
    qcrow = -1;
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
  int get_col() const {
    return ccol;
  }
  
  //Returns the row number of the QCore
  int get_row() const {
    return qcrow;
  }
  
  std::vector<bool> getHitmap();
    
  std::vector<bool> encodeQCore(bool is_new_col);
  
  const bool operator<(const QCore& other) {
    if (ccol==other.ccol) {
      return ccol < other.ccol;
    } else {
      return qcrow < other.qcrow;
    }
  }


 private:
  std::vector<bool> toRocCoordinates(std::vector<bool>& hitmap);

  std::vector<bool> int_to_binary(int num, int length);

  bool contains_hit(std::vector<bool>& hitmap);

  std::vector<bool> get_hitmap_code(std::vector<bool> hitmap);

};
#endif // QCORE_H
