#ifndef TRUTH_COUNTING_POLICY_H
#define TRUTH_COUNTING_POLICY_H

#include <cmath>
#include <cstdlib>

namespace TruthCountingPolicy
{
inline bool IsCountedChargedForActivity(long long pdg)
{
   const long long apdg = (pdg >= 0 ? pdg : -pdg);
   return (apdg == 11 || apdg == 13 || apdg == 15 ||
           apdg == 211 || apdg == 321 || apdg == 2212 ||
           apdg == 3112 || apdg == 3222 || apdg == 3312 ||
           apdg == 3334 || apdg == 411 || apdg == 431 ||
           apdg == 521 || apdg == 541 || apdg == 24);
}

inline double CountedChargeFromPdg(long long pdg)
{
   switch (pdg)
   {
      case 11:    return -1.0;
      case -11:   return +1.0;
      case 13:    return -1.0;
      case -13:   return +1.0;
      case 15:    return -1.0;
      case -15:   return +1.0;
      case 211:   return +1.0;
      case -211:  return -1.0;
      case 321:   return +1.0;
      case -321:  return -1.0;
      case 2212:  return +1.0;
      case -2212: return -1.0;
      case 3112:  return -1.0;
      case -3112: return +1.0;
      case 3222:  return +1.0;
      case -3222: return -1.0;
      case 3312:  return -1.0;
      case -3312: return +1.0;
      case 3334:  return -1.0;
      case -3334: return +1.0;
      case 411:   return +1.0;
      case -411:  return -1.0;
      case 431:   return +1.0;
      case -431:  return -1.0;
      case 521:   return +1.0;
      case -521:  return -1.0;
      case 541:   return +1.0;
      case -541:  return -1.0;
      case 24:    return +1.0;
      case -24:   return -1.0;
      default:    return 0.0;
   }
}

inline bool PassPIDFiducialFromMom(double px, double py, double pz,
                                   bool usePIDFiducial = true,
                                   double absCosMin = 0.15,
                                   double absCosMax = 0.675)
{
   if (!usePIDFiducial)
      return true;
   const double p2 = px * px + py * py + pz * pz;
   if (p2 <= 0.0)
      return false;
   const double absCos = std::abs(pz / std::sqrt(p2));
   return (absCos > absCosMin && absCos < absCosMax);
}

inline bool PassPtWindow(double px, double py,
                         double ptMin = 0.4,
                         double ptMax = 5.0)
{
   const double pt = std::sqrt(px * px + py * py);
   return (pt >= ptMin && pt < ptMax);
}

inline bool IsCountedPionForRatio(long long pdg, double px, double py, double pz,
                                  bool usePIDFiducial = true,
                                  double absCosMin = 0.15,
                                  double absCosMax = 0.675,
                                  double ptMin = 0.4,
                                  double ptMax = 5.0)
{
   return (std::abs(pdg) == 211 &&
           PassPtWindow(px, py, ptMin, ptMax) &&
           PassPIDFiducialFromMom(px, py, pz, usePIDFiducial, absCosMin, absCosMax));
}

inline bool IsCountedKaonForRatio(long long pdg, double px, double py, double pz,
                                  bool usePIDFiducial = true,
                                  double absCosMin = 0.15,
                                  double absCosMax = 0.675,
                                  double ptMin = 0.4,
                                  double ptMax = 5.0)
{
   return (std::abs(pdg) == 321 &&
           PassPtWindow(px, py, ptMin, ptMax) &&
           PassPIDFiducialFromMom(px, py, pz, usePIDFiducial, absCosMin, absCosMax));
}

inline bool IsCountedProtonForRatio(long long pdg, double px, double py, double pz,
                                    bool usePIDFiducial = true,
                                    double absCosMin = 0.15,
                                    double absCosMax = 0.675,
                                    double ptMin = 0.4,
                                    double ptMax = 5.0)
{
   return (std::abs(pdg) == 2212 &&
           PassPtWindow(px, py, ptMin, ptMax) &&
           PassPIDFiducialFromMom(px, py, pz, usePIDFiducial, absCosMin, absCosMax));
}
}

#endif
