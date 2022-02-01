#ifndef HIT_H
#define HIT_H

class Hit {
	private:

		int row_;
		int col_;
		int adc_;

	public:
		Hit(int row_num, int col_num, int adc_num);

		int row() const {

			return row_;
		}

		int col() const {
			return col_;
		}

		int adc() const {
			return adc_;
		}
};

#endif
