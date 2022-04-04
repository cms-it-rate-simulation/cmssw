#ifndef DataFormats_Phase2TrackerDigi_Hit_h
#define DataFormats_Phase2TrackerDigi_Hit_h

class Hit {
	public:
		Hit(int row, int col, int adc);

		int row() const {
			return row_;
		}

		int col() const {
			return col_;
		}

		int adc() const {
			return adc_;
		}

		void addADC(int adc);

	private:
		int row_;
		int col_;
		int adc_;
};

#endif
