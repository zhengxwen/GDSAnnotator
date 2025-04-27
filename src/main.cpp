// ===========================================================
//
// main.cpp: the C++ codes for the GDSAnnotator package
//
// Copyright (C) 2025    Xiuwen Zheng
//
// This file is part of GDSAnnotator.
//
// GDSAnnotator is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// GDSAnnotator is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GDSAnnotator.
// If not, see <http://www.gnu.org/licenses/>.

#include <R_GDS_CPP.h>
#include <Rinternals.h>
#include <string>
#include <vector>
#include <algorithm>


using namespace std;

// ===========================================================
// Library Functions
// ===========================================================

extern "C"
{

// Binary search function
// Returns the index of the target if found, otherwise -1
inline static int BinarySearch(const int arr[], int size, int target)
{
	int left = 0;
	int right = size - 1;
	while (left <= right)
	{
		int mid = left + (right - left) / 2;
		if (arr[mid] == target)
		{
			return mid;  // Target found
		} else if (arr[mid] < target)
		{
			left = mid + 1;  // Search in the right half
		} else {
			right = mid - 1;  // Search in the left half
		}
	}
	return -1;  // Target not found
}

LibExport SEXP SEQ_Find_Position(SEXP ptr, SEXP node_allele, SEXP pos,
	SEXP allele)
{
COREARRAY_TRY

	// initialize
	// positions in the gds file
	const int *gds_pos = (const int*)R_ExternalPtrAddr(ptr);
	const int n_gds_pos = Rf_asInteger(R_ExternalPtrProtected(ptr));
	// positions in the target
	const int *t_pos = INTEGER(pos);
	const int n_t_pos = Rf_length(pos);
	// the output vector
	rv_ans = PROTECT(NEW_INTEGER(n_t_pos));
	int *p_out_idx = INTEGER(rv_ans);
	//
	PdGDSObj nd = GDS_R_SEXP2Obj(node_allele, TRUE);

	vector<string> cur_allele;
	int cur_idx_st=0, cur_idx_ed=0;
	for (int p=gds_pos[cur_idx_st]; cur_idx_ed < n_gds_pos; cur_idx_ed++)
	{
		if (gds_pos[cur_idx_ed] != p) break;
	}

	// for-loop
	for (int i=0; i < n_t_pos; i++)
	{
		p_out_idx[i] = NA_INTEGER;
		if (t_pos[i] >= 0)
		{
			if (gds_pos[cur_idx_st] != t_pos[i])
			{
				// find the position via binary search
				int j = BinarySearch(gds_pos + cur_idx_ed,
					n_gds_pos - cur_idx_ed, t_pos[i]);
				if (j >= 0)
				{
					cur_idx_st = cur_idx_ed + j;
					const int p = gds_pos[cur_idx_st];
					// find the bound
					for (; cur_idx_st >= 0; cur_idx_st--)
						if (gds_pos[cur_idx_st] != p) break;
					cur_idx_st ++;
					cur_idx_ed += j;
					for (; cur_idx_ed < n_gds_pos; cur_idx_ed++)
						if (gds_pos[cur_idx_ed] != p) break;
					// need to reload cur_allele
					cur_allele.clear();
				} else
					continue;
			}
			if (cur_allele.empty())
			{
				// read alleles based on cur_idx_st & cur_idx_ed
				C_Int32 Start[1] = { cur_idx_st };
				C_Int32 Length[1] = { cur_idx_ed - cur_idx_st };
				cur_allele.resize(Length[0]);
				GDS_Array_ReadData(nd, Start, Length, &cur_allele[0],
					svStrUTF8);
			}
			// find the matched alleles
			string a = CHAR(STRING_ELT(allele, i));
			const vector<string>::iterator it =
				std::find(cur_allele.begin(), cur_allele.end(), a);
			if (it != cur_allele.end())
			{
				// 1-based
				p_out_idx[i] = it - cur_allele.begin() + cur_idx_st + 1;
			}
		}
	}

	UNPROTECT(1);

COREARRAY_CATCH
}


/// Initialize the package
LibExport void R_init_GDSAnnotator(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }

	static R_CallMethodDef callMethods[] =
	{
		CALL(SEQ_Find_Position, 4),
		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	Init_GDS_Routines();
}

} // extern "C"
