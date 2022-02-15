#ifndef HIT_H
#define HIT_H

class Hit {
	private:
		int row_;
		int col_;
		int adc_;

	public:
		Hit(int row_num, int col_num, int adc_num);
		void set_row(int row) { row_ = row; }
		void set_col(int col) { col_ = col; }
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
