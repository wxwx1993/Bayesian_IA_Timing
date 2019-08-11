// Generated by rstantools.  Do not edit by hand.

/*
    optimIA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    optimIA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with optimIA.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_commensurate_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_commensurate");
    reader.add_event(47, 45, "end", "model_commensurate");
    return reader;
}
#include <stan_meta_header.hpp>
class model_commensurate : public prob_grad {
private:
        int N;
        vector_d Y;
        std::vector<int> k;
        std::vector<int> hist;
public:
    model_commensurate(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_commensurate(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_commensurate_namespace::model_commensurate";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 1);
            current_statement_begin__ = 4;
            validate_non_negative_index("Y", "N", N);
            context__.validate_dims("data initialization", "Y", "vector_d", context__.to_vec(N));
            Y = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                Y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("k", "N", N);
            context__.validate_dims("data initialization", "k", "int", context__.to_vec(N));
            k = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("k");
            pos__ = 0;
            size_t k_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < k_k_0_max__; ++k_0__) {
                k[k_0__] = vals_i__[pos__++];
            }
            size_t k_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < k_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "k[i_0__]", k[i_0__], 1);
                check_less_or_equal(function__, "k[i_0__]", k[i_0__], 2);
            }
            current_statement_begin__ = 6;
            validate_non_negative_index("hist", "N", N);
            context__.validate_dims("data initialization", "hist", "int", context__.to_vec(N));
            hist = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("hist");
            pos__ = 0;
            size_t hist_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < hist_k_0_max__; ++k_0__) {
                hist[k_0__] = vals_i__[pos__++];
            }
            size_t hist_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < hist_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "hist[i_0__]", hist[i_0__], 0);
                check_less_or_equal(function__, "hist[i_0__]", hist[i_0__], 1);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 9;
            num_params_r__ += 1;
            current_statement_begin__ = 10;
            num_params_r__ += 1;
            current_statement_begin__ = 11;
            num_params_r__ += 1;
            current_statement_begin__ = 12;
            num_params_r__ += 1;
            current_statement_begin__ = 13;
            num_params_r__ += 1;
            current_statement_begin__ = 14;
            num_params_r__ += 1;
            current_statement_begin__ = 15;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_commensurate() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 9;
        if (!(context__.contains_r("mu01")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu01 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu01");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu01", "double", context__.to_vec());
        double mu01(0);
        mu01 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu01);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu01: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 10;
        if (!(context__.contains_r("mu02")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu02 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu02");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu02", "double", context__.to_vec());
        double mu02(0);
        mu02 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu02);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu02: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 11;
        if (!(context__.contains_r("s2")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable s2 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("s2");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "s2", "double", context__.to_vec());
        double s2(0);
        s2 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, s2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable s2: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 12;
        if (!(context__.contains_r("s3")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable s3 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("s3");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "s3", "double", context__.to_vec());
        double s3(0);
        s3 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, s3);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable s3: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 13;
        if (!(context__.contains_r("mu1")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu1 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu1");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu1", "double", context__.to_vec());
        double mu1(0);
        mu1 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu1: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 14;
        if (!(context__.contains_r("mu2")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu2 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu2");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu2", "double", context__.to_vec());
        double mu2(0);
        mu2 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu2: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 15;
        if (!(context__.contains_r("ss")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable ss missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("ss");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "ss", "double", context__.to_vec());
        double ss(0);
        ss = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, ss);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable ss: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 9;
            local_scalar_t__ mu01;
            (void) mu01;  // dummy to suppress unused var warning
            if (jacobian__)
                mu01 = in__.scalar_constrain(lp__);
            else
                mu01 = in__.scalar_constrain();
            current_statement_begin__ = 10;
            local_scalar_t__ mu02;
            (void) mu02;  // dummy to suppress unused var warning
            if (jacobian__)
                mu02 = in__.scalar_constrain(lp__);
            else
                mu02 = in__.scalar_constrain();
            current_statement_begin__ = 11;
            local_scalar_t__ s2;
            (void) s2;  // dummy to suppress unused var warning
            if (jacobian__)
                s2 = in__.scalar_lb_constrain(0, lp__);
            else
                s2 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 12;
            local_scalar_t__ s3;
            (void) s3;  // dummy to suppress unused var warning
            if (jacobian__)
                s3 = in__.scalar_lb_constrain(0, lp__);
            else
                s3 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 13;
            local_scalar_t__ mu1;
            (void) mu1;  // dummy to suppress unused var warning
            if (jacobian__)
                mu1 = in__.scalar_constrain(lp__);
            else
                mu1 = in__.scalar_constrain();
            current_statement_begin__ = 14;
            local_scalar_t__ mu2;
            (void) mu2;  // dummy to suppress unused var warning
            if (jacobian__)
                mu2 = in__.scalar_constrain(lp__);
            else
                mu2 = in__.scalar_constrain();
            current_statement_begin__ = 15;
            local_scalar_t__ ss;
            (void) ss;  // dummy to suppress unused var warning
            if (jacobian__)
                ss = in__.scalar_lb_constrain(0, lp__);
            else
                ss = in__.scalar_lb_constrain(0);
            // transformed parameters
            current_statement_begin__ = 18;
            validate_non_negative_index("mu", "4", 4);
            std::vector<local_scalar_t__> mu(4, local_scalar_t__(0));
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 19;
            local_scalar_t__ s2_sqrt;
            (void) s2_sqrt;  // dummy to suppress unused var warning
            stan::math::initialize(s2_sqrt, DUMMY_VAR__);
            stan::math::fill(s2_sqrt, DUMMY_VAR__);
            current_statement_begin__ = 20;
            local_scalar_t__ s3_sqrt;
            (void) s3_sqrt;  // dummy to suppress unused var warning
            stan::math::initialize(s3_sqrt, DUMMY_VAR__);
            stan::math::fill(s3_sqrt, DUMMY_VAR__);
            current_statement_begin__ = 21;
            local_scalar_t__ ss_sqrt;
            (void) ss_sqrt;  // dummy to suppress unused var warning
            stan::math::initialize(ss_sqrt, DUMMY_VAR__);
            stan::math::fill(ss_sqrt, DUMMY_VAR__);
            current_statement_begin__ = 22;
            local_scalar_t__ mudiff;
            (void) mudiff;  // dummy to suppress unused var warning
            stan::math::initialize(mudiff, DUMMY_VAR__);
            stan::math::fill(mudiff, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 23;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        mu1, 
                        "assigning variable mu");
            current_statement_begin__ = 24;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        mu2, 
                        "assigning variable mu");
            current_statement_begin__ = 25;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        mu01, 
                        "assigning variable mu");
            current_statement_begin__ = 26;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(4), stan::model::nil_index_list()), 
                        mu02, 
                        "assigning variable mu");
            current_statement_begin__ = 27;
            stan::math::assign(s2_sqrt, stan::math::sqrt(s2));
            current_statement_begin__ = 28;
            stan::math::assign(s3_sqrt, stan::math::sqrt(s3));
            current_statement_begin__ = 29;
            stan::math::assign(ss_sqrt, stan::math::sqrt(ss));
            current_statement_begin__ = 30;
            stan::math::assign(mudiff, (get_base1(mu, 2, "mu", 1) - get_base1(mu, 1, "mu", 1)));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 18;
            size_t mu_k_0_max__ = 4;
            for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
                if (stan::math::is_uninitialized(mu[k_0__])) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: mu" << "[" << k_0__ << "]";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 19;
            if (stan::math::is_uninitialized(s2_sqrt)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: s2_sqrt";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable s2_sqrt: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            current_statement_begin__ = 20;
            if (stan::math::is_uninitialized(s3_sqrt)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: s3_sqrt";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable s3_sqrt: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            current_statement_begin__ = 21;
            if (stan::math::is_uninitialized(ss_sqrt)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: ss_sqrt";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable ss_sqrt: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            current_statement_begin__ = 22;
            if (stan::math::is_uninitialized(mudiff)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: mudiff";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mudiff: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            // model body
            current_statement_begin__ = 34;
            lp_accum__.add(normal_log<propto__>(get_base1(mu, 3, "mu", 1), 0, 100));
            current_statement_begin__ = 35;
            lp_accum__.add(normal_log<propto__>(get_base1(mu, 4, "mu", 1), 0, 100));
            current_statement_begin__ = 36;
            lp_accum__.add(inv_gamma_log<propto__>(s2, 0.02, 1));
            current_statement_begin__ = 37;
            lp_accum__.add(inv_gamma_log<propto__>(s3, 0.02, 1));
            current_statement_begin__ = 38;
            lp_accum__.add(inv_gamma_log<propto__>(ss, 0.01, 1));
            current_statement_begin__ = 39;
            lp_accum__.add(normal_log<propto__>(get_base1(mu, 1, "mu", 1), get_base1(mu, 3, "mu", 1), s2_sqrt));
            current_statement_begin__ = 40;
            lp_accum__.add(normal_log<propto__>(get_base1(mu, 2, "mu", 1), get_base1(mu, 4, "mu", 1), s3_sqrt));
            current_statement_begin__ = 42;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 44;
                lp_accum__.add(normal_log<propto__>(get_base1(Y, i, "Y", 1), get_base1(mu, ((get_base1(hist, i, "hist", 1) * 2) + get_base1(k, i, "k", 1)), "mu", 1), ss_sqrt));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("mu01");
        names__.push_back("mu02");
        names__.push_back("s2");
        names__.push_back("s3");
        names__.push_back("mu1");
        names__.push_back("mu2");
        names__.push_back("ss");
        names__.push_back("mu");
        names__.push_back("s2_sqrt");
        names__.push_back("s3_sqrt");
        names__.push_back("ss_sqrt");
        names__.push_back("mudiff");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(4);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_commensurate_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double mu01 = in__.scalar_constrain();
        vars__.push_back(mu01);
        double mu02 = in__.scalar_constrain();
        vars__.push_back(mu02);
        double s2 = in__.scalar_lb_constrain(0);
        vars__.push_back(s2);
        double s3 = in__.scalar_lb_constrain(0);
        vars__.push_back(s3);
        double mu1 = in__.scalar_constrain();
        vars__.push_back(mu1);
        double mu2 = in__.scalar_constrain();
        vars__.push_back(mu2);
        double ss = in__.scalar_lb_constrain(0);
        vars__.push_back(ss);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 18;
            validate_non_negative_index("mu", "4", 4);
            std::vector<double> mu(4, double(0));
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 19;
            double s2_sqrt;
            (void) s2_sqrt;  // dummy to suppress unused var warning
            stan::math::initialize(s2_sqrt, DUMMY_VAR__);
            stan::math::fill(s2_sqrt, DUMMY_VAR__);
            current_statement_begin__ = 20;
            double s3_sqrt;
            (void) s3_sqrt;  // dummy to suppress unused var warning
            stan::math::initialize(s3_sqrt, DUMMY_VAR__);
            stan::math::fill(s3_sqrt, DUMMY_VAR__);
            current_statement_begin__ = 21;
            double ss_sqrt;
            (void) ss_sqrt;  // dummy to suppress unused var warning
            stan::math::initialize(ss_sqrt, DUMMY_VAR__);
            stan::math::fill(ss_sqrt, DUMMY_VAR__);
            current_statement_begin__ = 22;
            double mudiff;
            (void) mudiff;  // dummy to suppress unused var warning
            stan::math::initialize(mudiff, DUMMY_VAR__);
            stan::math::fill(mudiff, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 23;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        mu1, 
                        "assigning variable mu");
            current_statement_begin__ = 24;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        mu2, 
                        "assigning variable mu");
            current_statement_begin__ = 25;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        mu01, 
                        "assigning variable mu");
            current_statement_begin__ = 26;
            stan::model::assign(mu, 
                        stan::model::cons_list(stan::model::index_uni(4), stan::model::nil_index_list()), 
                        mu02, 
                        "assigning variable mu");
            current_statement_begin__ = 27;
            stan::math::assign(s2_sqrt, stan::math::sqrt(s2));
            current_statement_begin__ = 28;
            stan::math::assign(s3_sqrt, stan::math::sqrt(s3));
            current_statement_begin__ = 29;
            stan::math::assign(ss_sqrt, stan::math::sqrt(ss));
            current_statement_begin__ = 30;
            stan::math::assign(mudiff, (get_base1(mu, 2, "mu", 1) - get_base1(mu, 1, "mu", 1)));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t mu_k_0_max__ = 4;
                for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
                    vars__.push_back(mu[k_0__]);
                }
                vars__.push_back(s2_sqrt);
                vars__.push_back(s3_sqrt);
                vars__.push_back(ss_sqrt);
                vars__.push_back(mudiff);
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_commensurate";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu01";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu02";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "s2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "s3";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu1";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "ss";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t mu_k_0_max__ = 4;
            for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            param_name_stream__.str(std::string());
            param_name_stream__ << "s2_sqrt";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "s3_sqrt";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "ss_sqrt";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "mudiff";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu01";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu02";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "s2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "s3";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu1";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "ss";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t mu_k_0_max__ = 4;
            for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            param_name_stream__.str(std::string());
            param_name_stream__ << "s2_sqrt";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "s3_sqrt";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "ss_sqrt";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "mudiff";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_commensurate_namespace::model_commensurate stan_model;
#endif
