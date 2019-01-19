// DECLARATION OF stereoInt interfaces for STEREO
void stereoInt(Fimage * in1, Fimage * in2,
	       Fimage * disp, Fimage * cost, Fimage * dispr,
	       Fimage * costr, int patch_sz[], int hor_range[],
	       int ver_range[], char *method);
void stereoIntSlice(Fimage * in1, Fimage * in2, Fimage * disp,
		    Fimage * cost, Fimage * dispr, Fimage * costr,
		    int patch_sz[], int hor_range[], int ver_range[],
		    char *method);
void stereoIntSliceSubpix(Fimage * in1, Fimage * in2, Fimage * disp,
			  Fimage * cost, Fimage * dispr, Fimage * costr,
			  int patch_sz[], int hor_range[], int ver_range[],
			  char *method);

// DECLARATION OF bmIntSlice FOR 2D BLOCK MATCHING
void bmIntSlice(Fimage * in1, Fimage * in2,
		Fimage * disp, Fimage * cost, Fimage * dispr,
		Fimage * costr, int patch_sz[], int hor_range[],
		int ver_range[], char *method);

void bmIntSliceSubpix(Fimage * in1, Fimage * in2,
		      Fimage * disp, Fimage * cost, Fimage * dispr,
		      Fimage * costr, int patch_sz[], int hor_range[],
		      int ver_range[], char *method);
