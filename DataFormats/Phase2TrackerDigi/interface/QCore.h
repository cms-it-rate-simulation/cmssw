#ifndef QCORE_H
#define QCORE_H
#include<vector>

class QCore{

    public:
        QCore(
              int ccol_in,
              int qcrow_in,
              bool isneighbour_in,
              bool islast_in,
              std::vector<int> adcs_in
              );
        std::vector<int> adcs;
        bool islast;
        bool isneighbour;
        int ccol;
        int qcrow;

	int get_col();

	int get_row();

	std::vector<bool> to_ROC_coordinates(std::vector<bool> hitmap);
 
	std::vector<bool> get_hitmap();

	std::vector<bool> int_to_binary(int num, int length);

	bool contains_hit(std::vector<bool> hitmap);

	std::vector<bool> get_hitmap_code(std::vector<bool> hitmap);

	std::vector<bool> encode_qcore(bool is_new_col);
};
#endif // QCORE_H
