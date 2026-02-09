/*
 * Common_Functions.cpp
 *
 *  Created on: Aug 21, 2020
 *      Author: rodney
 */

#include "Common_Functions.hpp"
#include <cctype>
#include <cstring>
#include <cstdint>

namespace {

	// Check for NaN/Inf via IEEE 754 bit pattern.
	// std::isfinite() is unreliable under -ffast-math (-ffinite-math-only).
	bool is_non_finite(double value)
	{
		uint64_t bits;
		memcpy(&bits, &value, sizeof(bits));
		return (bits & UINT64_C(0x7FF0000000000000)) == UINT64_C(0x7FF0000000000000);
	}

	double interpolate_or_extrapolate(double x, double x_left, double y_left, double x_right, double y_right)
	{
		if (x_right == x_left)
		{
			return y_left;
		}

		// Use log-space interpolation for positive values (typical for density/temperature profiles).
		if (y_left > 0.0 && y_right > 0.0)
		{
			double log_y_left = log10(y_left);
			double log_y_right = log10(y_right);
			double dlogydx = (log_y_right - log_y_left) / (x_right - x_left);
			return pow(10.0, log_y_left + dlogydx * (x - x_left));
		}

		double dydx = (y_right - y_left) / (x_right - x_left);
		return y_left + dydx * (x - x_left);
	}

	void validate_x_column(const string &filename, const vector<double> &x_values)
	{
		if (x_values.size() < 2)
		{
			cout << "ERROR: \"" << filename << "\" has fewer than two data rows.\n";
			exit(1);
		}

		for (size_t i=0; i<x_values.size(); i++)
		{
			if (is_non_finite(x_values[i]))
			{
				cout << "ERROR: Non-finite value in first column of \"" << filename << "\" at row " << (i+1) << ".\n";
				exit(1);
			}
			if (i > 0 && x_values[i] <= x_values[i-1])
			{
				cout << "ERROR: First column in \"" << filename << "\" must be strictly increasing.\n";
				exit(1);
			}
		}
	}

	void sanitize_non_finite_column(const string &filename, int col_index, const vector<double> &x_values, vector<double> &y_values)
	{
		if (y_values.size() != x_values.size())
		{
			cout << "ERROR: Mismatched column lengths while reading \"" << filename << "\" (column " << col_index << ").\n";
			exit(1);
		}

		vector<size_t> finite_indices;
		finite_indices.reserve(y_values.size());
		for (size_t i=0; i<y_values.size(); i++)
		{
			if (!is_non_finite(y_values[i]))
			{
				finite_indices.push_back(i);
			}
		}

		if (finite_indices.empty())
		{
			cout << "ERROR: Column " << col_index << " in \"" << filename << "\" contains no finite values.\n";
			exit(1);
		}

		int replaced = 0;

		// Fill leading non-finite values.
		size_t first_finite = finite_indices.front();
		if (first_finite > 0)
		{
			size_t right0 = finite_indices[0];
			size_t right1 = (finite_indices.size() > 1) ? finite_indices[1] : finite_indices[0];
			for (size_t i=0; i<first_finite; i++)
			{
				if (right0 == right1)
				{
					y_values[i] = y_values[right0];
				}
				else
				{
					y_values[i] = interpolate_or_extrapolate(x_values[i], x_values[right0], y_values[right0], x_values[right1], y_values[right1]);
				}
				replaced++;
			}
		}

		// Fill interior non-finite values.
		for (size_t k=0; k + 1 < finite_indices.size(); k++)
		{
			size_t left = finite_indices[k];
			size_t right = finite_indices[k+1];
			for (size_t i=left+1; i<right; i++)
			{
				if (is_non_finite(y_values[i]))
				{
					y_values[i] = interpolate_or_extrapolate(x_values[i], x_values[left], y_values[left], x_values[right], y_values[right]);
					replaced++;
				}
			}
		}

		// Fill trailing non-finite values by extrapolating from the last two valid points.
		// If only one finite value exists, fall back to clamping.
		size_t last_finite = finite_indices.back();
		if (last_finite + 1 < y_values.size())
		{
			size_t left = (finite_indices.size() > 1) ? finite_indices[finite_indices.size()-2] : last_finite;
			size_t right = last_finite;
			for (size_t i=last_finite+1; i<y_values.size(); i++)
			{
				if (left == right)
				{
					y_values[i] = y_values[right];
				}
				else
				{
					y_values[i] = interpolate_or_extrapolate(x_values[i], x_values[left], y_values[left], x_values[right], y_values[right]);
				}
				replaced++;
			}
		}

		// Final guard in case any non-finite values remain.
		for (size_t i=0; i<y_values.size(); i++)
		{
			if (is_non_finite(y_values[i]))
			{
				size_t left = i;
				while (left > 0 && is_non_finite(y_values[left])) left--;
				size_t right = i;
				while (right + 1 < y_values.size() && is_non_finite(y_values[right])) right++;

				if (!is_non_finite(y_values[left]) && !is_non_finite(y_values[right]) && right != left)
				{
					y_values[i] = interpolate_or_extrapolate(x_values[i], x_values[left], y_values[left], x_values[right], y_values[right]);
				}
				else if (!is_non_finite(y_values[left]))
				{
					y_values[i] = y_values[left];
				}
				else if (!is_non_finite(y_values[right]))
				{
					y_values[i] = y_values[right];
				}
				else
				{
					cout << "ERROR: Unable to sanitize non-finite value in \"" << filename << "\" column " << col_index << ".\n";
					exit(1);
				}
				replaced++;
			}
		}

		if (replaced > 0)
		{
			cout << "WARNING: Replaced " << replaced << " non-finite value(s) in \"" << filename << "\" column " << col_index << ".\n";
		}
	}

	void import_csv_generic(const string &filename, vector<vector<double>*> columns)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			cout << "\"" << filename << "\" not found!\n";
			exit(1);
		}

		size_t line_number = 0;
		string line;
		while (getline(infile, line))
		{
			line_number++;

			if (line.empty() || line[0] == '#' || std::all_of(line.begin(), line.end(), [](unsigned char c){ return std::isspace(c); }))
			{
				continue;
			}

			replace(line.begin(), line.end(), ',', ' ');
			stringstream str(line);

			// Peek at first token — if it isn't numeric, treat line as a header and skip.
			string first_token;
			if (!(str >> first_token)) continue;
			{
				string lower = first_token;
				for (auto &c : lower) c = std::tolower(c);
				size_t pos = 0;
				bool is_number = false;
				try { stod(first_token, &pos); is_number = (pos > 0); } catch (...) {}
				if (!is_number && lower != "nan" && lower != "-nan" && lower != "inf" && lower != "-inf")
				{
					continue;  // skip non-numeric header line
				}
			}
			// Reset stream to re-read the full line including the first token.
			str.clear();
			str.str(line);

			for (size_t i=0; i<columns.size(); i++)
			{
				string token;
				if (!(str >> token))
				{
					cout << "ERROR: Could not parse column " << (i+1) << " in \"" << filename << "\" at line " << line_number << ".\n";
					exit(1);
				}

				// Handle nan/NaN/NAN and inf explicitly.
				// numeric_limits<double>::quiet_NaN() is unreliable under -ffast-math
				// (clang++ can optimise it to garbage), so build NaN/Inf from bit patterns.
				string lower_token = token;
				for (auto &c : lower_token) c = std::tolower(c);
				double v;
				if (lower_token == "nan" || lower_token == "-nan")
				{
					uint64_t nan_bits = UINT64_C(0x7FF8000000000000);
					memcpy(&v, &nan_bits, sizeof(v));
				}
				else if (lower_token == "inf" || lower_token == "infinity")
				{
					uint64_t inf_bits = UINT64_C(0x7FF0000000000000);
					memcpy(&v, &inf_bits, sizeof(v));
				}
				else if (lower_token == "-inf" || lower_token == "-infinity")
				{
					uint64_t ninf_bits = UINT64_C(0xFFF0000000000000);
					memcpy(&v, &ninf_bits, sizeof(v));
				}
				else
				{
					try { v = stod(token); }
					catch (...)
					{
						cout << "ERROR: Could not parse column " << (i+1) << " in \"" << filename << "\" at line " << line_number << ".\n";
						exit(1);
					}
				}
				columns[i]->push_back(v);
			}
		}
		infile.close();

		if (columns.empty() || columns[0]->empty())
		{
			cout << "ERROR: \"" << filename << "\" contains no readable data rows.\n";
			exit(1);
		}

		validate_x_column(filename, *(columns[0]));
		for (size_t i=1; i<columns.size(); i++)
		{
			sanitize_non_finite_column(filename, static_cast<int>(i+1), *(columns[0]), *(columns[i]));
		}
	}
}

// function to check if custom random seed exists in local file "rng_seed"
// if file does not exist, uses system clock to generate seed
static long long get_seed()
{
	long long s = 0;
	ifstream seedfile;
	seedfile.open("rng_seed");
	if (!seedfile.good())
	{
		s = chrono::high_resolution_clock::now().time_since_epoch().count();
	}
	else
	{
		string line;
		getline(seedfile, line);
		s = stoll(line);
	}
	seedfile.close();
	return s;
}

// seed random number generator using get_seed() function above
// to access externally, must include "Common_Functions.hpp" and call
// using common::get_rand() (will return uniform real between 0 and 1)
static long long seed = get_seed();
static mt19937 rand_generator(seed);   // Mersenne Twister PRNG (apparently, pretty good)
static uniform_real_distribution<double> rand_dist(0.0, 1.0);  // dist to be used with get_rand()

namespace constants {
	const double pi    = M_PI;            // pi [unitless]
	const double twopi = 2*pi;             // 2*pi [unitless]
	const double k_b   = 1.380649e-16;     // Boltzmann's Constant [erg/K]
	const double c     = 29979245800.0;    // Speed of Light in Vacuum [cm/s]
	const double G     = 6.67430e-8;       // Gravitational Constant [cm^3/g/s^2]
	const double amu   = 1.660538782e-24;  // Atomic Mass Unit [g]
	const double m_e   = 9.10938215e-28;   // Electron Mass [g]
	const double q_e   = 1.602176487e-19;  // Elementary Charge [C]
	const double jev   = q_e;              // Joules/Electron Volt [unitless]
	const double ergev = jev*1.0e7;        // ergs/Electron Volt [unitless]
}

namespace common {

	void import_csv(string filename, vector<double> &col1, vector<double> &col2)
	{
		import_csv_generic(filename, {&col1, &col2});
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3)
	{
		import_csv_generic(filename, {&col1, &col2, &col3});
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4)
	{
		import_csv_generic(filename, {&col1, &col2, &col3, &col4});
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5)
	{
		import_csv_generic(filename, {&col1, &col2, &col3, &col4, &col5});
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5, vector<double> &col6)
	{
		import_csv_generic(filename, {&col1, &col2, &col3, &col4, &col5, &col6});
	}

	// returns interpolated value at x from parallel arrays (x_data, y_data)
	// assumes that x_data has at least two elements, is sorted and is strictly monotonic increasing
	double interpolate(vector<double> &x_data, vector<double> &y_data, double x)
	{
		/*
		int size = x_data.size();

		int i = 0;                    // find left end of interval for interpolation
		if (x >= x_data[size - 2])    // special case: beyond right end
		{
			i = size - 2;
		}
		else
		{
			while (x > x_data[i+1]) i++;
		}
		double x_left = x_data[i], y_left = y_data[i], x_right = x_data[i+1], y_right = y_data[i+1];  // points on either side (unless beyond ends)
		*/

		if (x <= x_data[0])
		{
			return y_data[0];
		}
		else if (x >= x_data.back())
		{
			return y_data.back();
		}
		else
		{
			int i = lower_bound(x_data.begin(), x_data.end(), x) - x_data.begin();
			double x_left = x_data[i-1], x_right = x_data[i], y_left = y_data[i-1], y_right = y_data[i];
			double dydx = (y_right - y_left) / (x_right - x_left);   // gradient
			return y_left + dydx * (x - x_left);   // linear interpolation
		}
	}

	// returns interpolated value at x from arrays x_data, y_data, when y_data is log-scaled
	double interpolate_logy(vector<double> &x_data, vector<double> &y_data, double x)
	{
		if (x <= x_data[0])
		{
			return y_data[0];
		}
		else if (x >= x_data.back())
		{
			return y_data.back();
		}
		else
		{
			int i = lower_bound(x_data.begin(), x_data.end(), x) - x_data.begin();
			double x_left = x_data[i-1], x_right = x_data[i], y_left = log10(y_data[i-1]), y_right = log10(y_data[i]);
			double dydx = (y_right - y_left) / (x_right - x_left);   // gradient
			return pow(10.0, y_left + dydx * (x - x_left));   // linear interpolation
		}
	}

	// returns uniformly distributed random number from interval [0, 1)
	double get_rand()
	{
		return rand_dist(rand_generator);
	}

	// returns uniformly distributed random integer between lower and upper (inclusive)
	int get_rand_int(int lower, int upper)
	{
		uniform_int_distribution<int> dist(lower, upper);
		return dist(rand_generator);
	}
}
