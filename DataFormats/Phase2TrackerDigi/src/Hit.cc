#include "../interface/Hit.h"

//Describes a hit (meant for use in 4x4 sensor coordinates) by its row number, column number, and adc value
Hit::Hit(int row_num, int col_num, int adc_num) {
	row_ = row_num;
	col_ = col_num;
	adc_ = adc_num;
}

//Adds the input value to the adc to a max of 15 - used to combine two hits that are at the same position
void Hit::addADC(int adc) {
	adc_ += adc;

	if(adc_ > 15) {
		adc_ = 15;
	}
}
