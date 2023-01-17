// MathLibrary.cpp : Defines the exported functions for the DLL.
#include "pch.h" // use pch.h in Visual Studio 2019 and later
#include "dllmain.h"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <array>
#include <numeric>
#include <math.h>
#include <limits>

#include <objscip/objscip.h>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

// TODO to be removed later
#include <scip/struct_cons.h>
#include <scip/type_cons.h>
#include <scip/cons_linear.h>

using namespace std;

double const INVALID_VALUE = 0.0;

/* returns the ratio */
class Ratio {
public:
    double comp(double number_);
    Ratio(double base_);
    Ratio& setBase(double newbase_);
private:
    double base;
};
Ratio::Ratio(double base_) : base(base_) {
    //   initialization is enough :)
}

//compute the ratio between number and base
double Ratio::comp(double number_) {
    return number_ / max(base, 1.0);
}

class VecStats {
public:
    double mean(void);
    double mmin(void);
    double mmax(void);
    double std(void);
    double median(void);
    double density(void);
    VecStats(vector<double>& vec_, SCIP* scip_);

private:
    vector<double>& vec;
    SCIP* scip;
};

VecStats::VecStats(vector<double>& vec_, SCIP* scip_) : vec(vec_), scip(scip_) {

}

/**
   Calculates maximal element of a vector.

   @return Maximal element of the vector or INVALID_VALUE if vector is empty.
*/
double VecStats::mmax() {
    if (vec.size() == 0) return INVALID_VALUE;
    return *std::max_element(begin(vec), end(vec));
}

/**
   Calculates minimal element of a vector.

   @return Minimal element of the vector or INVALID_VALUE if vector is empty.
*/
double VecStats::mmin() {
    if (vec.size() == 0) return INVALID_VALUE;
    return *std::min_element(begin(vec), end(vec));
}

/**
   Calculates mean element of a vector.

   @return Mean value of the vector or INVALID_VALUE if vector is empty.
*/
double VecStats::mean() {
    double vecsum = 0.0;
    if (vec.size() == 0) return INVALID_VALUE;
    for (size_t j = 0; j < vec.size(); ++j)
        vecsum += vec[j];

    return vecsum / (double)vec.size();
}

/**
   Calculates median element of a vector.

   @return Median element of the vector or INVALID_VALUE if vector is empty.
*/
double VecStats::std() {
    //   numerically stable computation of variance
    size_t n;
    double mean = 0.0;
    double M2 = 0.0;
    if (vec.size() == 0) return INVALID_VALUE;

    /* iterate over vector elements */
    for (n = 0; n < vec.size(); ++n)
    {
        double delta = vec[n] - mean;
        mean += delta / ((double)n + 1.0);
        double delta2 = vec[n] - mean;
        M2 += delta * delta2;
    }

    if (n < 2)
        return 0.0;
    else
        return M2 / ((double)n - 1.0);
}

/**
   Calculates median element of a vector.

   @return Median element of the vector or INVALID_VALUE if vector is empty.
*/
double VecStats::median() {
    if (vec.empty()) return INVALID_VALUE;

    vector<double> veccopy(vec);

    /* set element on index vector_size / 2 to element that would appear at
    that position after full sorting of vector elements */
    nth_element(veccopy.begin(), veccopy.begin() + veccopy.size() / 2, veccopy.end());
    double median = veccopy[veccopy.size() / 2];

    /* for vectors with even number of elements, set element on index vector_size / 2 - 1
    to element that would appear there after vector sorting */
    if (veccopy.size() % 2 == 0) {
        nth_element(veccopy.begin(), veccopy.begin() + veccopy.size() / 2 - 1, veccopy.end());
        median += veccopy[veccopy.size() / 2 - 1];
        median /= 2;
    }

    return median;
}

/**
   Calculates vector density.

   @return Percent of nonzero elements in vector.
*/
double VecStats::density() {
    if (vec.empty()) return INVALID_VALUE;

    return count_if(vec.begin(), vec.end(), [&](double x) { return !SCIPisZero(scip, x); }) / (1.0 * vec.size());
}

Ratio& Ratio::setBase(double newbase_) {
    base = newbase_;
    return *this;
}
const char* delim(void) {
    return ",";
}

/**
   Calculates sum of all vector elements.

   @param   vec   Reference to the vector whose sum will be calculated.

   @return Vector sum.
*/
int vecsum(std::vector<int>& vec) {
    int sum = 0;
    /* loop over vector and compute the sum */
    for (int j = 0; j < vec.size(); ++j)
        sum += vec[j];
    return sum;
}

/**
   Prints vector statistics including mean, minimum, maximum, standard deviation, median and vector density.

   @param   vec            Reference to a vector for which statistics is calculated.
   @param   scip           Pointer to the instance of SCIP solver.
   @param   print_density  Flag that control if vector density will be printed or not (if true it will be printed)
   @param   print_last_delim Flag to control printing of last delimiter (if false, no delimiter will end the last element)
*/
void writeVecStats(std::ofstream& features, const char* feature_name, std::vector<double>& vec, SCIP* scip, bool print_density = false) {
    VecStats stats{ vec, scip };

    features << "mean_" << feature_name << delim() << stats.mean() << endl;
    features << "min_" << feature_name << delim() << stats.mmin() << endl;
    features << "max_" << feature_name << delim() << stats.mmax() << endl;
    features << "std_" << feature_name << delim() << stats.std() << endl;
    features << "median_" << feature_name << delim() << stats.median() << endl;
    if (print_density) features << "densiry_" << feature_name << delim() << stats.density() << endl;
}


/**
 * Normalizes constraint coefficients, including handsides.
 *
 *   @param   rhs               Reference to constraint right hand side.
 *   @param   lhs               Reference to constraint left hand side.
 *   @param   constraintvals    Vector of constraint coefficients.
 */
void normalize(SCIP_Real& rhs, SCIP_Real& lhs, vector<SCIP_Real>& constraintvals, SCIP_Real infinity) {
    if (constraintvals.empty()) return;

    vector<SCIP_Real> absconstraintvals;
    transform(constraintvals.begin(), constraintvals.end(), back_inserter(absconstraintvals), [](double x) { return abs(x); });
    double normalization_factor = *max_element(absconstraintvals.begin(), absconstraintvals.end());
    if (normalization_factor != 0) {
        /* normalize variable coefficients */
        for_each(constraintvals.begin(), constraintvals.end(), [=](double& x) { x /= normalization_factor; });

        /* normalize hand sides */
        if (rhs < infinity)
            rhs /= normalization_factor;
        if (lhs > -infinity)
            lhs /= normalization_factor;
    }
}

enum VARTYPE { VARTYPE_FIXED, VARTYPE_BINARY, VARTYPE_INTEGER, VARTYPE_IMPLINTEGER, VARTYPE_CONTINUOUS };

/**
   Determine type of objective function variable.

   @param   scip        Pointer to an instance of SCIP solver.
   @param   scip_var    Pointer to objective function variable.

   @return Type of scip_var (type from VARTYPE enumeration)
*/
VARTYPE varGetType(SCIP* scip, SCIP_VAR* scip_var) {
    if (SCIPisEQ(scip, SCIPvarGetLbGlobal(scip_var), SCIPvarGetUbGlobal(scip_var)) == TRUE) {
        return VARTYPE_FIXED;
    }
    else {
        auto scip_vartype = SCIPvarGetType(scip_var);
        switch (scip_vartype)
        {
        case SCIP_VARTYPE_BINARY:     return VARTYPE_BINARY;
        case SCIP_VARTYPE_INTEGER:    return VARTYPE_INTEGER;
        case SCIP_VARTYPE_IMPLINT:    return VARTYPE_IMPLINTEGER;
        case SCIP_VARTYPE_CONTINUOUS: return VARTYPE_CONTINUOUS;
        default: return VARTYPE_CONTINUOUS;
        }
    }
}

/**
   Calculates number of fixed variables
   where fixed variables are variables with equal upper an lower bounds.

   @param   scip  Pointer to an instance of SCIP solver.
   @return  Number of fixed variables
*/
int getNFixedVars(SCIP* scip) {
    int nvars = SCIPgetNVars(scip);
    SCIP_VAR** scip_vars = SCIPgetVars(scip);

    return count_if(scip_vars, scip_vars + nvars, [&](SCIP_VAR* scip_var) {
        return (SCIPisEQ(scip, SCIPvarGetLbGlobal(scip_var), SCIPvarGetUbGlobal(scip_var)) == TRUE);
        });
}



/* update clique entry vector for the current constraint */
void updateCliqueInformation(
    SCIP* scip,
    SCIP_CONS* cons,
    std::vector<double>& cliquentries
)
{
    SCIP_Real* vals = SCIPgetValsLinear(scip, cons);
    int nvals = SCIPgetNVarsLinear(scip, cons);
    SCIP_VAR** vars = SCIPgetVarsLinear(scip, cons);
    SCIP_Real rhs = SCIPgetRhsLinear(scip, cons);
    SCIP_Real lhs = SCIPgetLhsLinear(scip, cons);
    int i;
    SCIP_Real a_min = SCIPinfinity(scip);
    SCIP_Real a_negsum = 0;
    SCIP_Real c_min = 0;
    SCIP_Real a_min2 = SCIPinfinity(scip);
    int nbin = 0;

    /* in the case of an equation or a ranged row, treat both sides as individual cliques */
    while (!SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs))
    {
        SCIP_Bool negate = FALSE;

        /* negate constraints of the form a^T x >= b as -a^T x <= -b */
        if (SCIPisInfinity(scip, rhs))
        {
            rhs = -lhs;
            lhs = -SCIPinfinity(scip);
            negate = TRUE;
        }
        /* get
         *  - a_min the absolute minimum binary coefficient
         *  - a_negsum the sum of negative binary coefficients
         *  - c_min the minimum activity contribution of all non-binary variables
         *  - the second smallest binary coefficient a_min2
         *  - nbin the number of binary variables
         */
        double* zero = 0;

        for (i = 0; i < nvals; ++i)
        {
            SCIP_Real val = negate ? -vals[i] : vals[i];

            if (SCIPisZero(scip, val))
                continue;

            if (!SCIPvarIsBinary(vars[i]))
            {
                SCIP_Real minbound = vals > zero ? SCIPvarGetLbGlobal(vars[i]) : SCIPvarGetUbGlobal(vars[i]);

                if (SCIPisInfinity(scip, REALABS(minbound)))
                    c_min = -SCIPinfinity(scip);
                else if (!SCIPisInfinity(scip, -c_min))
                    c_min += val * minbound;
            }
            else
            {
                ++nbin;
                a_negsum += (val < 0 ? val : 0.0);

                /* update minimum and second smallest absolute binary coefficients */
                if (a_min > REALABS(val))
                {
                    a_min2 = a_min;
                    a_min = REALABS(val);
                }
                else if (a_min2 > REALABS(val))
                    a_min2 = REALABS(val);
            }
        }

        /* if d - c_min - a_negsum - a_min < a_min2, the binary variables of this constraint form a clique */
        if (nbin > 0 && (SCIPisLT(scip, rhs - c_min - a_negsum - a_min, a_min2)))
        {
            /* add the binary density of this clique to the cliqueentries */
            cliquentries.push_back(nbin / (double)SCIPgetNBinVars(scip));
        }

        /* relax right hand side */
        rhs = SCIPinfinity(scip);
    }
}

int write_features(const char* problem_filename, const char* csv_filename, const char* instance_name) {
    SCIP* scip;

    SCIP_CALL(SCIPcreate(&scip));
    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 0));
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));

    /*SCIPsetMessagehdlrQuiet(scip, TRUE);*/
    SCIPreadProb(scip, problem_filename, NULL);

    cout.precision(numeric_limits<double>::max_digits10);
    // disable upgrading to more specific constraints because of clique function
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/upgrade/indicator", FALSE));
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/upgrade/knapsack", FALSE));
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/upgrade/logicor", FALSE));
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/upgrade/setppc", FALSE));
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/upgrade/varbound", FALSE));
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/upgrade/xor", FALSE));

    SCIPpresolve(scip);

    /******************************************************/
    /*             Open csv-file for writeing             */
    /******************************************************/
    std::ofstream features;
    features.open(csv_filename);

    /******************************************************/
    /*             Features calculation                   */
    /******************************************************/

    features << "INSTANCE_NAME" << delim() << instance_name << endl; /* INSTANCE_NAME,instance_name */
    features << "LOG_VARS" << delim() << SCIPgetNVars(scip) << endl; /* LOG_VARS,log_vars */

    /* binary, integer, continuous, implicit integer and fixed variables */
    features << "BIN_VARS" << delim() << SCIPgetNBinVars(scip) << endl; /* BIN_VARS,bin_var */
    features << "INT_VARS" << delim() << SCIPgetNIntVars(scip) << endl; /* INT_VARS,int_var */
    features << "CONT_VARS" << delim() << SCIPgetNContVars(scip) << endl; /* CONT_VARS,cont_var */
    features << "IMPLINT_VARS" << delim() << SCIPgetNImplVars(scip) << endl; /* IMPLINT_VARS, implint_vars */
    features << "FIXED_VARS" << delim() << getNFixedVars(scip) << endl; /* FIXED_VARS, fixed_vars */

    /* collect objective statistics */
    std::vector<int> objsByType(5);
    std::vector<double> absobjnonzerovals;
    std::vector<double> upperbounds;
    std::vector<double> lowerbounds;
    std::vector<double> boundranges;
    std::array< std::vector<double>, 6> objsarr;
    for (int j = 0; j < SCIPgetNVars(scip); ++j)
    {
        SCIP_VAR** scip_vars = SCIPgetVars(scip);
        double obj = SCIPvarGetObj(scip_vars[j]);
        int increment = !SCIPisZero(scip, obj);
        objsarr[0].push_back(obj);

        /* required for objective dynamism */
        if (!SCIPisZero(scip, obj))
            absobjnonzerovals.push_back(abs(obj));

        /* switch through the types */
        switch (varGetType(scip, scip_vars[j]))
        {
        case VARTYPE_BINARY:
            objsarr[1].push_back(obj);
            objsByType[0] += increment;
            break;
        case VARTYPE_INTEGER:
            objsarr[2].push_back(obj);
            objsByType[1] += increment;
            break;
        case VARTYPE_CONTINUOUS:
            objsarr[3].push_back(obj);
            objsByType[2] += increment;
            break;
        case VARTYPE_IMPLINTEGER:
            objsarr[4].push_back(obj);
            objsByType[3] += increment;
            break;
        case VARTYPE_FIXED:
            objsarr[5].push_back(obj);
            objsByType[4] += increment;
            break;
        default:
            break;
        }

        SCIP_Real upperbound = SCIPvarGetUbGlobal(scip_vars[j]);
        SCIP_Real lowerbound = SCIPvarGetLbGlobal(scip_vars[j]);
        if (!SCIPisInfinity(scip, REALABS(upperbound)))
            upperbounds.push_back(upperbound);
        if (!SCIPisInfinity(scip, REALABS(lowerbound)))
            lowerbounds.push_back(lowerbound);

        if (!SCIPisInfinity(scip, REALABS(upperbound)) && !SCIPisInfinity(scip, REALABS(lowerbound)))
            boundranges.push_back(upperbound - lowerbound);

    }

    /* number of variables (all, typed) with nonzero objective coefficients */
    features << "TOTAL_OBJNONZERO_VARS" << delim() << vecsum(objsByType) << endl; /* TOTAL_OBJNONZERO_VARS,tot_objnonzero */
    features << "BIN_OBJNONZERO_VARS" << delim() << objsByType[0] << endl; /* BIN_OBJNONZERO_VARS,bin_objnonzero */
    features << "INT_OBJNONZERO_VARS" << delim() << objsByType[1] << endl; /* INT_OBJNONZERO_VARS,int_objnonzero */
    features << "CONT_OBJNONZERO_VARS" << delim() << objsByType[2] << endl; /* CONT_OBJNONZERO_VARS,cont_objnonzero */
    features << "IMPLINT_OBJNONZERO_VARS" << delim() << objsByType[3] << endl; /* IMPLINT_OBJNONZERO_VARS, implint_objnonzero */
    features << "FIXED_OBJNONZERO_VARS" << delim() << objsByType[4] << endl; /* FIXED_OBJNONZERO_VARS, fixed_objnonzero */

    /* author bzfgojic: variable lower/upper bound statistics */
    Ratio ratio(SCIPgetNVars(scip));
    features << "PER_FINITE_UB" << delim() << ratio.comp(upperbounds.size()) << endl; /* PER_FINITE_UB,per_finit_upper */
    features << "PER_FINITE_LB" << delim() << ratio.comp(lowerbounds.size()) << endl; /* PER_FINITE_LB,per_finit_lower */
    writeVecStats(features, "UPPERBOUNDS", upperbounds, scip, true);   // also prints vector density when second parameter is true /* UPPERBOUNDS_*,ub_* */
    writeVecStats(features, "LOWERBOUNDS", lowerbounds, scip, true);   /* LOWERBOUNDS_*,lb_* */
    writeVecStats(features, "BOUNDRANGE", boundranges, scip, false);  /* BOUNDRANGE_*, boundrange_* density is unnecessary */

    /* calculating min and max element of absolute nonzero objective coefficients */
    if (!absobjnonzerovals.empty()) {
        auto minmax = minmax_element(absobjnonzerovals.begin(), absobjnonzerovals.end());
        assert(!SCIPisZero(scip, *minmax.first));
        features << "OBJECTIVE_DYNAMISM" << delim() << log10((*minmax.second) / (*minmax.first)) << endl;
    }
    else {
        features << "OBJECTIVE_DYNAMISM" << delim() << INVALID_VALUE << endl;
    }

    /* objective vector|{all, bin, int, cont, implint, fixed} : mean, min, max, std, median */
    for (int j = 0; j < 6; ++j)
    {
        writeVecStats(features, "OBJCOEFF", objsarr[j], scip); /* OBJCOEFF_*,objcoeff_* */
    }

    /* number of constraints */
    SCIP_CONS** conss = SCIPgetConss(scip);
    int nconss = SCIPgetNConss(scip);
    int nindicatorconss = 0;
    int nsosconss = 0;
    int notherconss = 0;
    int nrangedrows = 0;

    features << "LOG_CONSTR" << delim() << nconss << endl; /* LOG_CONSTR,log_constr */

    double nnonzeroconsscoef = 0;

    SCIP_CONSHDLR* linear = SCIPfindConshdlr(scip, "linear");
    SCIP_CONSHDLR* indicator = SCIPfindConshdlr(scip, "indicator");
    SCIP_CONSHDLR* sos = SCIPfindConshdlr(scip, "SOS1");
    assert(sos != nullptr);
    assert(linear != nullptr);
    assert(indicator != nullptr);

    std::vector<SCIP_Real> sides;
    std::vector<SCIP_Real> righthandsides;
    std::vector<SCIP_Real> lefthandsides;
    std::vector<SCIP_Real> consnvars;
    std::vector<SCIP_Real> consmincoefs;
    std::vector<SCIP_Real> consmaxcoefs;
    std::vector<SCIP_Real> consmeancoefs;
    std::vector<SCIP_Real> consstdcoefs;
    std::vector<SCIP_Real> consrangecoefs;
    std::vector<SCIP_Real> conslogmaxminratio;      // log max/min absolute value of nonzero coeffs for all constraints
    std::vector<double> cliquentries;

    for (int i = 0; i < nconss; ++i)
    {
        SCIP_CONS* cons = conss[i];
        SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);

        /* most of the constraint classification only works for linear constraints */
        if (conshdlr != linear)
        {
            if (conshdlr == sos)
                ++nsosconss;
            else if (conshdlr == indicator)
                ++nindicatorconss;
            else
                ++notherconss;

            continue;
        }


        if (conshdlr == linear)
        {

            SCIP_Real rhs = SCIPgetRhsLinear(scip, cons);
            SCIP_Real lhs = SCIPgetLhsLinear(scip, cons);

            /* update the clique information based on this constraint */
            updateCliqueInformation(scip, cons, cliquentries);

            int nconsvars = SCIPgetNVarsLinear(scip, cons);             // number of variables included in the constraint
            vector<SCIP_Real> consvals;
            if (nconsvars) {
                SCIP_Real* consvalsarr = SCIPgetValsLinear(scip, cons);  // values of the coeff in the  linear constraint
                consvals.assign(consvalsarr, consvalsarr + nconsvars);
                normalize(rhs, lhs, consvals, SCIPinfinity(scip));
            }

            if (!SCIPisInfinity(scip, rhs))
            {
                sides.push_back(std::abs(rhs));
                righthandsides.push_back(rhs);
            }
            if (!SCIPisInfinity(scip, -lhs))
            {
                lefthandsides.push_back(lhs);
                if (!SCIPisEQ(scip, lhs, rhs)) {
                    sides.push_back(std::abs(lhs));

                    if (!SCIPisInfinity(scip, rhs))
                        nrangedrows++;
                }
            }

            if (nconsvars) {
                nnonzeroconsscoef += count_if(consvals.begin(), consvals.end(), [&](double x) { return !SCIPisZero(scip, x); }); // number of all nonzero coeff for current constraint
                consnvars.push_back(nconsvars);

                vector<SCIP_Real> absconsvals;
                std::transform(consvals.begin(), consvals.end(), std::back_inserter(absconsvals), [](double x) { return std::abs(x); });

                /* author bzfgojic: calculating log max/min, max required for dynamism */
                SCIP_Real consmaxabscoef, consminabscoef;
                auto new_end = remove_if(absconsvals.begin(), absconsvals.end(), [&](double x) { return SCIPisZero(scip, x); });
                absconsvals.resize(new_end - absconsvals.begin());

                auto minmax = minmax_element(absconsvals.begin(), absconsvals.end());
                if (SCIPisPositive(scip, *minmax.first)) {
                    conslogmaxminratio.push_back(log10((*minmax.second) / (*minmax.first)));
                }

                VecStats consstats{ consvals, scip };

                double min = consstats.mmin();
                double max = consstats.mmax();

                consmeancoefs.push_back(consstats.mean()); // vector of mean values of the abs values of the coeff 
                consmincoefs.push_back(min);
                consmaxcoefs.push_back(max);
                consstdcoefs.push_back(consstats.std());
                consrangecoefs.push_back(Ratio{ max }.comp(min));
            }
        }
    }

    /* linear constraint classification */
    SCIP_LINCONSSTATS linconsstats;
    SCIP_LINCONSSTATS* linconsstatsptr;
    linconsstatsptr = &linconsstats;

    SCIP_CALL(SCIPclassifyConstraintTypesLinear(scip, linconsstatsptr));

    int nlinempty = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_EMPTY);
    int nlinfree = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_FREE);
    int nlinsing = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_SINGLETON);
    int nlinaggr = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_AGGREGATION);
    int nlinprec = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_PRECEDENCE);
    int nlinvarbd = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_VARBOUND);
    int nlinsetpart = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_SETPARTITION);
    int nlinsetpack = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_SETPACKING);
    int nlinsetcov = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_SETCOVERING);
    int nlincard = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_CARDINALITY);
    int nlininvknap = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_INVKNAPSACK);
    int nlineqknap = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_EQKNAPSACK);
    int nlinbinpack = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_BINPACKING);
    int nlinknap = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_KNAPSACK);
    int nlinintknap = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_INTKNAPSACK);
    int nlinmixbin = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_MIXEDBINARY);
    int nlingen = SCIPlinConsStatsGetTypeCount(linconsstatsptr, SCIP_LINCONSTYPE_GENERAL);

    features << "LINTOTAL_CONSTR" << delim() << SCIPlinConsStatsGetSum(linconsstatsptr) << endl;
    features << "LINEMPTY_CONSTR" << delim() << nlinempty << endl; /* LINEMPTY_CONSTR, linempty_constr */
    features << "LINFREE_CONSTR" << delim() << nlinfree << endl;
    features << "LINSING_CONSTR" << delim() << nlinsing << endl;
    features << "LINAGGR_CONSTR" << delim() << nlinaggr << endl;
    features << "LINPREC_CONSTR" << delim() << nlinprec << endl;
    features << "LINVARBD_CONSTR" << delim() << nlinvarbd << endl;
    features << "LINSETPART_CONSTR" << delim() << nlinsetpart << endl;
    features << "LINSETPACK_CONSTR" << delim() << nlinsetpack << endl;
    features << "LINSETPACK_CONSTR" << delim() << nlinsetpack << endl;
    features << "LINSETCOV_CONSTR" << delim() << nlinsetcov << endl;
    features << "LINCARD_CONSTR" << delim() << nlincard << endl;
    features << "LININVKNAP_CONSTR" << delim() << nlininvknap << endl;
    features << "LINEQKNAP_CONSTR" << delim() << nlineqknap << endl;
    features << "LINBINPACK_CONSTR" << delim() << nlinbinpack << endl;
    features << "LINBINPACK_CONSTR" << delim() << nlinknap << endl;
    features << "LINKNAPSACK_CONSTR" << delim() << nlinknap << endl;
    features << "LININTKNAPSACK_CONSTR" << delim() << nlinintknap << endl;
    features << "LINMIXBIN_CONSTR" << delim() << nlinmixbin << endl;
    features << "LINGEN_CONSTR" << delim() << nlingen << endl;

    writeVecStats(features, "CONSTR", sides, scip); /* CONSTR_*,constr_* */

    /* author bzfgojic: statistics for constraint lefthand and righthandsides vectors */
    writeVecStats(features, "RH_CONSTR", righthandsides, scip, true); /* RH_CONSTR_*,rh_constr_* */
    writeVecStats(features, "LH_CONSTR", lefthandsides, scip, true); /* LH_CONSTR_*,lh_constr_* */

    features << "RANGED_ROWS" << delim() << nrangedrows << endl; /* RANGED_ROWS, ranged_rows */

    /* author bzfgojic: ratio of finite constraint sides (right and left) and total number of specific sides */
    ratio.setBase(nconss);
    features << "RH_CONSTR_RATIO" << delim() << ratio.comp(righthandsides.size()) << endl;
    features << "LH_CONSTR_RATIO" << delim() << ratio.comp(lefthandsides.size()) << endl;

    /* author bzfhende: variable number per constraint |{all, bin, int, cont} mean, min, max, std */
    writeVecStats(features, "consnvars", consnvars, scip); /* NVAR_ALL_*,nvar_all_* */

    /* author bzfhende: coefficient features */
    writeVecStats(features, "CONSTR_COEFF_MEAN", consmeancoefs, scip); /* CONSTR_COEFF_MEAN_*,constr_coeff_mean_* */
    writeVecStats(features, "CONSTR_COEFF_MIN", consmincoefs, scip); /* CONSTR_COEFF_MIN_*,constr_coeff_min_* */
    writeVecStats(features, "CONSTR_COEFF_MAX", consmaxcoefs, scip); /* CONSTR_COEFF_MAX_*,constr_coeff_max_* */
    writeVecStats(features, "CONSTR_COEFF_STD", consstdcoefs, scip); /* CONSTR_COEFF_STD_*,constr_coeff_std_* */
    writeVecStats(features, "CONSTR_COEFF_RATIO", consrangecoefs, scip); /* CONSTR_COEFF_RATIO_*,constr_coeff_ratio_* */

    /* author bzfgojic: constraint matrix density */
    double base = (double)nconss;
    base *= (double)SCIPgetNVars(scip);
    ratio.setBase(base);
    features << "CONSTR_MATRIX_DENSITY" << delim() << ratio.comp(nnonzeroconsscoef) << endl; /* CONSTR_MATRIX_DENSITY,constr_matrix_density */
    features << "NONZEROES" << delim() << nnonzeroconsscoef << endl; /* NONZEROES, nonzeroes */

    features << "INDICATOR_CONSTR" << delim() << nindicatorconss << endl; /* INDICATOR_CONSTR,indicator_constr */
    features << "SOS_CONSTR" << delim() << nsosconss << endl; /* SOS_CONSTR, sos_constr  */
    features << "OTHER_CONSTR" << delim() << notherconss << endl; /* OTHER_CONSTR,other_constr */

    /* clique statistics */
    writeVecStats(features, "CLIQUE", cliquentries, scip); /* CLIQUE_*,clique_* */

    /* dynamism statistics. Use false to indicate that no delimiter is necessary */
    writeVecStats(features, "DYNAMISM", conslogmaxminratio, scip, false); /* DYNAMISM_*,dynamism_* */

    /* author bzfhende: TODO constraint number per variable | {all, bin, int, cont} mean, min, max, std */
    SCIP_CALL(SCIPfree(&scip));
    features.close();
    return 0;
}
