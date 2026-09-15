#ifndef ROTATIONAL_CROSS_SECTIONS_HPP_
#define ROTATIONAL_CROSS_SECTIONS_HPP_

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Cross sections are tabulated at TOTAL energy K + E(ji), in eV.
// CSV columns: total_energy_eV,ji,jf,sigma_cm2 (including elastic diagonals).
struct RotationalTransition {
    int ji, jf;
    double sigma, delta_e;
};

class RotationalCrossSections {
public:
    struct Curve {
        int ji, jf;
        std::vector<double> energy, sigma, reduced;
        size_t source = 0;
    };
    std::vector<Curve> curves;
    double B = 0.0;
    int max_ji = 0;
    // Thermal sensitivity option: Boltzmann weights only for j <= this level,
    // renormalized over the kept levels (-1 keeps the full distribution).
    int population_max_j = -1;
    std::vector<bool> diagonals;

    void load(const std::string &path, double rot_B) {
        if (!(rot_B > 0.0)) throw std::runtime_error("State cross sections require rot_const_cm1 > 0");
        B = rot_B;
        max_ji = 0;
        std::ifstream input(path);
        if (!input) throw std::runtime_error("Cannot open state cross sections: " + path);
        std::map<std::pair<int,int>, Curve> data;
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0]=='#' || line.find("total_energy_eV")==0) continue;
            std::replace(line.begin(),line.end(),',',' ');
            std::istringstream row(line);
            double e,s; int i,j;
            if (!(row>>e>>i>>j>>s) || !std::isfinite(e) || !std::isfinite(s) ||
                e<0 || s<0 || i<0 || j<0)
                throw std::runtime_error("Invalid state cross-section row in " + path);
            auto &c=data[{i,j}]; c.ji=i; c.jf=j;
            if (!c.energy.empty() && e<=c.energy.back())
                throw std::runtime_error("State cross-section energies must increase per channel: " + path);
            c.energy.push_back(e); c.sigma.push_back(s); max_ji=std::max(max_ji,i);
        }
        curves.clear();
        for (auto &entry:data) curves.push_back(std::move(entry.second));
        if (curves.empty()) throw std::runtime_error("Empty state cross sections: " + path);
        if (data.find({0,0})==data.end()) throw std::runtime_error("State cross sections must include the ground elastic channel");
        diagonals.assign(max_ji+1,false);
        for (const auto &c:curves) if (c.ji==c.jf) diagonals[c.ji]=true;
        std::map<std::pair<int,int>,size_t> indices;
        for (size_t k=0;k<curves.size();++k) indices[{curves[k].ji,curves[k].jf}]=k;
        for (size_t k=0;k<curves.size();++k) {
            auto &c=curves[k]; c.source=k;
            const double ei=level(c.ji), ef=level(c.jf);
            for (size_t n=0;n<c.energy.size();++n) {
                double r=c.sigma[n];
                if (c.ji<c.jf) {
                    if (c.energy[n]<=ef) {
                        if (r>0) throw std::runtime_error("Positive cross section below channel opening");
                        r=0;
                    } else r*=(2*c.ji+1)*(c.energy[n]-ei)/(c.energy[n]-ef);
                } else if (c.ji>c.jf) r*=2*c.ji+1;
                c.reduced.push_back(r);
            }
            auto reverse=indices.find({c.jf,c.ji});
            if (c.ji>c.jf && reverse!=indices.end()) c.source=reverse->second;
        }
    }

    double level(int j) const { return B*double(j)*(j+1.0); }

    // Interpolate a shared reduced cross section R(E) for each reversible pair:
    // g_low*K_low*sigma_up = (E-E_high)*R = g_high*K_high*sigma_down.
    // This preserves detailed balance between nodes and during extrapolation.
    // Hold R outside the grid (elastic sigma is held). This is an explicit
    // threshold/extrapolation model, not a measured low-energy cross section.
    double cross_section(const Curve &c, double kinetic) const {
        const double ei=level(c.ji), ef=level(c.jf);
        if (!(kinetic>0.0) || kinetic<=ef-ei) return 0.0;
        const double e=kinetic+ei;
        const auto &source=curves[c.source];
        double r;
        if (e<=source.energy.front()) r=source.reduced.front();
        else if (e>=source.energy.back()) r=source.reduced.back();
        else {
            size_t hi=std::upper_bound(source.energy.begin(),source.energy.end(),e)-source.energy.begin();
            double t=(e-source.energy[hi-1])/(source.energy[hi]-source.energy[hi-1]);
            r=(1-t)*source.reduced[hi-1]+t*source.reduced[hi];
        }
        if (c.ji==c.jf) return r;
        return r/(2*c.ji+1)*(c.ji<c.jf?(kinetic-(ef-ei))/kinetic:1.0);
    }

    std::vector<double> populations(double temperature, bool thermal, const std::string &species) const {
        std::vector<double> p(max_ji+1,0.0);
        if (!thermal) { p[0]=1.0; return p; }
        if (max_ji==0) throw std::runtime_error("Thermal rotations require excited-initial-state cross sections for " + species);
        if (!(temperature>0.0)) throw std::runtime_error("Thermal rotations require positive temperature");
        const bool capped=population_max_j>=0;
        if (population_max_j>max_ji) throw std::runtime_error("rot_population_max_j exceeds tabulated initial states for " + species);
        double z=0,covered=0;
        for (int j=0; ; ++j) {
            double spin=1;
            if (species=="CO2" && j%2) spin=0; // 12C16O2, v=0
            if (species=="N2") spin=(j%2)?3:6; // 14N2
            double w=spin*(2.0*j+1)*std::exp(-level(j)/(8.617333262145e-5*temperature));
            z+=w;
            if (j<=max_ji && diagonals[j]) { p[j]=w; covered+=w; }
            else if (capped && w>0) throw std::runtime_error("Missing initial state below rot_population_max_j for " + species);
            if (capped && j==population_max_j) break;
            if (j>max_ji && level(j)>40*8.617333262145e-5*temperature) break;
            if (j>100000) throw std::runtime_error("Invalid rotational temperature");
        }
        if (!capped && covered/z<0.999) throw std::runtime_error("Rotational basis covers less than 99.9% of population for " + species);
        for (double &x:p) x/=covered;
        return p;
    }

    std::vector<RotationalTransition> evaluate(double kinetic, double temperature, bool thermal,
                                             const std::string &species, double ground_elastic) const {
        auto p=populations(temperature,thermal,species);
        std::vector<RotationalTransition> result;
        for (const auto &c:curves) {
            if (p[c.ji]==0) continue;
            double s=(c.ji==0 && c.jf==0)?ground_elastic:cross_section(c,kinetic);
            // Exact zero for elastic channels: FMA contraction can leave +-1e-19 eV in E(j)-E(j).
            if (s>0) result.push_back({c.ji,c.jf,p[c.ji]*s,c.ji==c.jf?0.0:level(c.jf)-level(c.ji)});
        }
        return result;
    }

    // Bound valid for every kinetic energy and any normalized population.
    double upper_bound(double ground_elastic_max) const {
        std::map<int,double> bound;
        for (const auto &c:curves) {
            const auto &source=curves[c.source];
            double value=*std::max_element(source.reduced.begin(),source.reduced.end());
            if (c.ji!=c.jf) value/=2*c.ji+1;
            bound[c.ji]+=(c.ji==0 && c.jf==0)?ground_elastic_max:value;
        }
        double largest=ground_elastic_max;
        for (const auto &b:bound) largest=std::max(largest,b.second);
        return largest*(1+1e-12);
    }
};
#endif
