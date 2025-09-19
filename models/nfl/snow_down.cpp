#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <snow_down.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  styr.allocate("styr");
  endyr.allocate("endyr");
  dat_yr.allocate("dat_yr");
  years.allocate(1,dat_yr,"years");
  size_n.allocate("size_n");
  sizes.allocate(1,size_n,"sizes");
  imm_n_obs.allocate(styr,endyr,"imm_n_obs");
  imm_n_size_obs.allocate(styr,endyr,1,size_n,"imm_n_size_obs");
  mat_n_obs.allocate(styr,endyr,"mat_n_obs");
  mat_n_size_obs.allocate(styr,endyr,1,size_n,"mat_n_size_obs");
  prop_term_molt.allocate(styr,endyr,1,size_n,"prop_term_molt");
  size_trans.allocate(1,size_n,1,size_n,"size_trans");
  prop_skip.allocate(styr,endyr,1,size_n,"prop_skip");
  sigma_numbers_imm.allocate("sigma_numbers_imm");
  sigma_numbers_mat.allocate("sigma_numbers_mat");
  mat_eff_samp.allocate("mat_eff_samp");
  imm_eff_samp.allocate("imm_eff_samp");
  log_mu_m_prior.allocate("log_mu_m_prior");
  log_mu_m_mat_prior.allocate(1,2,"log_mu_m_mat_prior");
  est_m_devs.allocate("est_m_devs");
  est_q_devs.allocate("est_q_devs");
  est_m_mat_devs.allocate("est_m_mat_devs");
  est_q_mat_devs.allocate("est_q_mat_devs");
  est_sigma_m.allocate("est_sigma_m");
  sigma_m_mu.allocate(1,3,"sigma_m_mu");
  smooth_q_weight.allocate("smooth_q_weight");
  smooth_m_weight.allocate(1,2,"smooth_m_weight");
  est_log_m_mu.allocate("est_log_m_mu");
  est_log_m_mat_mu.allocate("est_log_m_mat_mu");
  est_sigma_q.allocate("est_sigma_q");
  est_m_lg_devs.allocate("est_m_lg_devs");
  smooth_f_weight.allocate("smooth_f_weight");
  smooth_surv_weight.allocate("smooth_surv_weight");
cout<<"prop_skip"<<prop_skip<<endl; 
cout<<"sigma_numbers_imm"<<sigma_numbers_imm<<endl; 
 ad_comm::change_datafile_name("catch_dat.DAT");
  ret_cat_yr_n.allocate("ret_cat_yr_n");
  cat_st.allocate("cat_st");
  cat_end.allocate("cat_end");
  ret_cat_yrs.allocate(1,ret_cat_yr_n,"ret_cat_yrs");
  ret_cat_numbers.allocate(cat_st,cat_end,"ret_cat_numbers");
  disc_cat_yr_n.allocate("disc_cat_yr_n");
  d_cat_st.allocate("d_cat_st");
  d_cat_end.allocate("d_cat_end");
  disc_cat_yrs.allocate(1,disc_cat_yr_n,"disc_cat_yrs");
  disc_cat_numbers.allocate(d_cat_st,d_cat_end,"disc_cat_numbers");
  sigma_numbers_ret.allocate("sigma_numbers_ret");
  sigma_numbers_disc.allocate("sigma_numbers_disc");
  discard_survival.allocate("discard_survival");
  ret_eff_samp.allocate("ret_eff_samp");
  disc_eff_samp.allocate("disc_eff_samp");
  ret_cat_size.allocate(d_cat_st,d_cat_end,1,size_n,"ret_cat_size");
  disc_cat_size.allocate(d_cat_st,d_cat_end,1,size_n,"disc_cat_size");
cout<<"sigma_numbers_ret"<<sigma_numbers_ret<<endl;
cout<<"disc_eff_samp"<<disc_eff_samp<<endl;
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_n_imm.allocate(1,size_n,1.01,30,1,"log_n_imm");
  log_n_mat.allocate(1,size_n,1.01,30,1,"log_n_mat");
  nat_m_dev.allocate(styr,endyr,-4,4,est_m_devs,"nat_m_dev");
  nat_m_mat_dev.allocate(styr,endyr,-4,4,est_m_mat_devs,"nat_m_mat_dev");
  nat_m_lg_dev.allocate(styr,endyr,-4,4,est_m_mat_devs,"nat_m_lg_dev");
  q_dev.allocate(styr,endyr,-0.2,0.2,est_q_devs,"q_dev");
  q_mat_dev.allocate(styr,endyr,-0.2,0.2,est_q_mat_devs,"q_mat_dev");
  q_lg_dev.allocate(styr,endyr,-0.2,0.2,est_q_mat_devs,"q_lg_dev");
  log_avg_rec.allocate(1,40,"log_avg_rec");
  rec_devs.allocate(styr,endyr,-10,10,1,"rec_devs");
  sigma_m.allocate(1,3,0.01,4,est_sigma_m,"sigma_m");
  log_m_mu.allocate(1-5,3,est_log_m_mu,"log_m_mu");
  log_m_mat_mu.allocate(1,2,-5,3,est_log_m_mat_mu,"log_m_mat_mu");
  prop_rec.allocate(1,2,0.00001,200,"prop_rec");
  sigma_q.allocate(1,2,0.01,4,est_sigma_q,"sigma_q");
  log_f.allocate(-5,5,"log_f");
  f_dev.allocate(styr,endyr,-5,5,"f_dev");
  fish_ret_sel_50.allocate(25,150,3,"fish_ret_sel_50");
  fish_ret_sel_slope.allocate(0.0001,20,3,"fish_ret_sel_slope");
  fish_tot_sel_50.allocate(25,150,2,"fish_tot_sel_50");
  fish_tot_sel_slope.allocate(0.0001,20,2,"fish_tot_sel_slope");
  surv_sel_50.allocate(25,150,2,"surv_sel_50");
  surv_sel_slope.allocate(0.0001,20,2,"surv_sel_slope");
  imm_n_size_pred.allocate(styr,endyr,1,size_n,"imm_n_size_pred");
  #ifndef NO_AD_INITIALIZE
    imm_n_size_pred.initialize();
  #endif
  mat_n_size_pred.allocate(styr,endyr,1,size_n,"mat_n_size_pred");
  #ifndef NO_AD_INITIALIZE
    mat_n_size_pred.initialize();
  #endif
  nat_m.allocate(styr,endyr,1,size_n,"nat_m");
  #ifndef NO_AD_INITIALIZE
    nat_m.initialize();
  #endif
  nat_m_mat.allocate(styr,endyr,1,size_n,"nat_m_mat");
  #ifndef NO_AD_INITIALIZE
    nat_m_mat.initialize();
  #endif
  selectivity.allocate(styr,endyr,1,size_n,"selectivity");
  #ifndef NO_AD_INITIALIZE
    selectivity.initialize();
  #endif
  selectivity_mat.allocate(styr,endyr,1,size_n,"selectivity_mat");
  #ifndef NO_AD_INITIALIZE
    selectivity_mat.initialize();
  #endif
  surv_sel.allocate(1,size_n,"surv_sel");
  #ifndef NO_AD_INITIALIZE
    surv_sel.initialize();
  #endif
  total_fish_sel.allocate(1,size_n,"total_fish_sel");
  #ifndef NO_AD_INITIALIZE
    total_fish_sel.initialize();
  #endif
  retain_fish_sel.allocate(1,size_n,"retain_fish_sel");
  #ifndef NO_AD_INITIALIZE
    retain_fish_sel.initialize();
  #endif
  pred_retained_n.allocate(styr,endyr,"pred_retained_n");
  #ifndef NO_AD_INITIALIZE
    pred_retained_n.initialize();
  #endif
  pred_discard_n.allocate(styr,endyr,"pred_discard_n");
  #ifndef NO_AD_INITIALIZE
    pred_discard_n.initialize();
  #endif
  pred_retained_size_comp.allocate(styr,endyr-1,1,size_n,"pred_retained_size_comp");
  #ifndef NO_AD_INITIALIZE
    pred_retained_size_comp.initialize();
  #endif
  pred_discard_size_comp.allocate(styr,endyr-1,1,size_n,"pred_discard_size_comp");
  #ifndef NO_AD_INITIALIZE
    pred_discard_size_comp.initialize();
  #endif
  temp_imm.allocate(1,size_n,"temp_imm");
  #ifndef NO_AD_INITIALIZE
    temp_imm.initialize();
  #endif
  temp_skip.allocate(1,size_n,"temp_skip");
  #ifndef NO_AD_INITIALIZE
    temp_skip.initialize();
  #endif
  temp_noskip.allocate(1,size_n,"temp_noskip");
  #ifndef NO_AD_INITIALIZE
    temp_noskip.initialize();
  #endif
  temp_mat.allocate(1,size_n,"temp_mat");
  #ifndef NO_AD_INITIALIZE
    temp_mat.initialize();
  #endif
  trans_imm.allocate(1,size_n,"trans_imm");
  #ifndef NO_AD_INITIALIZE
    trans_imm.initialize();
  #endif
  temp_catch_imm.allocate(1,size_n,"temp_catch_imm");
  #ifndef NO_AD_INITIALIZE
    temp_catch_imm.initialize();
  #endif
  temp_catch_mat.allocate(1,size_n,"temp_catch_mat");
  #ifndef NO_AD_INITIALIZE
    temp_catch_mat.initialize();
  #endif
  sum_imm_numbers_obs.allocate(styr,endyr,"sum_imm_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_imm_numbers_obs.initialize();
  #endif
  sum_mat_numbers_obs.allocate(styr,endyr,"sum_mat_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_mat_numbers_obs.initialize();
  #endif
  imm_numbers_pred.allocate(styr,endyr,"imm_numbers_pred");
  #ifndef NO_AD_INITIALIZE
    imm_numbers_pred.initialize();
  #endif
  mat_numbers_pred.allocate(styr,endyr,"mat_numbers_pred");
  #ifndef NO_AD_INITIALIZE
    mat_numbers_pred.initialize();
  #endif
  sum_ret_numbers_obs.allocate(styr,endyr,"sum_ret_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_ret_numbers_obs.initialize();
  #endif
  sum_disc_numbers_obs.allocate(styr,endyr,"sum_disc_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_disc_numbers_obs.initialize();
  #endif
  imm_num_like.allocate("imm_num_like");
  #ifndef NO_AD_INITIALIZE
  imm_num_like.initialize();
  #endif
  mat_num_like.allocate("mat_num_like");
  #ifndef NO_AD_INITIALIZE
  mat_num_like.initialize();
  #endif
  ret_cat_like.allocate("ret_cat_like");
  #ifndef NO_AD_INITIALIZE
  ret_cat_like.initialize();
  #endif
  disc_cat_like.allocate("disc_cat_like");
  #ifndef NO_AD_INITIALIZE
  disc_cat_like.initialize();
  #endif
  imm_like.allocate("imm_like");
  #ifndef NO_AD_INITIALIZE
  imm_like.initialize();
  #endif
  mat_like.allocate("mat_like");
  #ifndef NO_AD_INITIALIZE
  mat_like.initialize();
  #endif
  ret_comp_like.allocate("ret_comp_like");
  #ifndef NO_AD_INITIALIZE
  ret_comp_like.initialize();
  #endif
  disc_comp_like.allocate("disc_comp_like");
  #ifndef NO_AD_INITIALIZE
  disc_comp_like.initialize();
  #endif
  use_term_molt.allocate(styr,endyr,1,size_n,"use_term_molt");
  #ifndef NO_AD_INITIALIZE
    use_term_molt.initialize();
  #endif
  nat_m_like.allocate("nat_m_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_like.initialize();
  #endif
  nat_m_mat_like.allocate("nat_m_mat_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mat_like.initialize();
  #endif
  nat_m_mu_like.allocate("nat_m_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mu_like.initialize();
  #endif
  nat_m_lg_like.allocate("nat_m_lg_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_lg_like.initialize();
  #endif
  nat_m_mat_mu_like.allocate("nat_m_mat_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mat_mu_like.initialize();
  #endif
  nat_m_lg_mu_like.allocate("nat_m_lg_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_lg_mu_like.initialize();
  #endif
  smooth_q_like.allocate("smooth_q_like");
  #ifndef NO_AD_INITIALIZE
  smooth_q_like.initialize();
  #endif
  smooth_m_like.allocate("smooth_m_like");
  #ifndef NO_AD_INITIALIZE
  smooth_m_like.initialize();
  #endif
  q_like.allocate("q_like");
  #ifndef NO_AD_INITIALIZE
  q_like.initialize();
  #endif
  q_mat_like.allocate("q_mat_like");
  #ifndef NO_AD_INITIALIZE
  q_mat_like.initialize();
  #endif
  smooth_f_like.allocate("smooth_f_like");
  #ifndef NO_AD_INITIALIZE
  smooth_f_like.initialize();
  #endif
  surv_sel_prior.allocate("surv_sel_prior");
  #ifndef NO_AD_INITIALIZE
  surv_sel_prior.initialize();
  #endif
  smooth_surv_like.allocate("smooth_surv_like");
  #ifndef NO_AD_INITIALIZE
  smooth_surv_like.initialize();
  #endif
  temp_prop_rec.allocate(1,3,"temp_prop_rec");
  #ifndef NO_AD_INITIALIZE
    temp_prop_rec.initialize();
  #endif
  tot_prop_rec.allocate("tot_prop_rec");
  #ifndef NO_AD_INITIALIZE
  tot_prop_rec.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
 dvariable fpen=0.0;
 for(int size=1;size<=size_n;size++)
  {
   imm_n_size_pred(styr,size) = exp(log_n_imm(size));
   mat_n_size_pred(styr,size) = exp(log_n_mat(size));
   total_fish_sel(size) = 1 / (1+exp(-fish_tot_sel_slope*(sizes(size)-fish_tot_sel_50))) ; 
   retain_fish_sel(size) = 1 / (1+exp(-fish_ret_sel_slope*(sizes(size)-fish_ret_sel_50)));
   surv_sel(size) = 1 / (1+exp(-surv_sel_slope*(sizes(size)-surv_sel_50)));
   }
 for(int year=styr;year<=endyr;year++)
  for(int size=1;size<=size_n;size++)
  {
  nat_m(year,size) = log_m_mu + nat_m_dev(year);
  nat_m_mat(year,size) = log_m_mat_mu(1) + nat_m_dev(year);
  selectivity(year,size) = surv_sel(size) + q_dev(year);
  selectivity_mat(year,size) = surv_sel(size) + q_dev(year);
  //==option for estimating separate devs by maturity state
  if(est_q_mat_devs>0)
   selectivity_mat(year,size) = surv_sel(size) + q_mat_dev(year);  
  if(est_m_mat_devs>0)
   nat_m_mat(year,size) = log_m_mat_mu(1) + nat_m_mat_dev(year);
  if(est_m_lg_devs>0 & size > 13)
  {
   nat_m(year,size) = log_m_mat_mu(2) + nat_m_lg_dev(year);
   nat_m_mat(year,size) = log_m_mat_mu(2) + nat_m_lg_dev(year);
  }
  }
  temp_prop_rec.initialize();
  tot_prop_rec = 20;
  for(int i=1;i<=2;i++)
   tot_prop_rec += prop_rec(i);
   temp_prop_rec(1) = 20/ tot_prop_rec;
  for(int i = 1;i<=2;i++)
   temp_prop_rec(i+1) = prop_rec(i)/tot_prop_rec;
 rec_devs(endyr) =0;
 // use_term_molt = prop_term_molt;
 // for(int size = 1; size<=size_n;size++)
 //  use_term_molt(2020,size) = prob_mat_2020(size);
  use_term_molt = prop_term_molt;
  //for(int size = 1; size<=16;size++)
  // use_term_molt(2020,size) = 1 / (1+exp(-prob_mat_2020(size)));
 //==make last year rec dev and f dev equal to the average for now
 // rec_devs(endyr) =0;
 // f_dev(endyr) = 0;
  // cout<<"log_n_imm"<<log_n_imm<<endl;
  // cout<<"log_n_mat"<<log_n_mat<<endl;
  // cout<<"nat_m_dev"<<nat_m_dev<<endl;
  // cout<<"nat_m_mat_dev"<<nat_m_mat_dev<<endl;
  // cout<<"log_avg_rec"<<log_avg_rec<<endl;
  // cout<<"rec_devs"<<rec_devs<<endl;
  // cout<<"log_m_mu"<<log_m_mu<<endl;
  // cout<<"prop_rec"<<prop_rec<<endl;
  // cout<<"log_f"<<log_f<<endl;
  // cout<<"f_dev"<<f_dev<<endl;
  // cout<<"fish_ret_sel_50"<<fish_ret_sel_50<<endl;
  // cout<<"fish_ret_sel_slope"<<fish_ret_sel_slope<<endl;
  // cout<<"fish_tot_sel_50"<<fish_tot_sel_50<<endl;
  // cout<<"fish_tot_sel_slope"<<fish_tot_sel_slope<<endl;
  // cout<<"surv_sel"<<surv_sel<<endl;
  // cout<<"temp_prop_rec"<<temp_prop_rec<<endl;
  // cout<<"nat_m"<<nat_m<<endl;
   // cout<<"selectivity"<<selectivity<<endl;
    // cout<<"retain_fish_sel"<<retain_fish_sel<<endl;
	 // cout<<"total_fish_sel"<<total_fish_sel<<endl;
    // cout<<"imm_n_size_pred"<<imm_n_size_pred(styr)<<endl;
	 // cout<<"mat_n_size_pred"<<mat_n_size_pred(styr)<<endl;	 
 for(int year=styr;year<endyr;year++)
  {
  for (int size=1;size<=size_n;size++) 
      {
		temp_imm(size) = imm_n_size_pred(year,size) * exp(-1*(0.59)*exp( nat_m(year,size)));
		temp_mat(size) = mat_n_size_pred(year,size) * exp(-1*(0.59)*exp( nat_m_mat(year,size)));
	   }	
	   for(int size=1;size<=size_n;size++)
	   {	   
	   temp_catch_imm(size) = temp_imm(size) * (1 -exp(-(exp(log_f+f_dev(year))*total_fish_sel(size))));
	   temp_catch_mat(size) = temp_mat(size) * (1 -exp(-(exp(log_f+f_dev(year))*total_fish_sel(size))));
	   pred_retained_size_comp(year,size) = temp_catch_imm(size)*retain_fish_sel(size) + temp_catch_mat(size)*retain_fish_sel(size); 
	   pred_discard_size_comp(year,size) = temp_catch_imm(size)*(1-retain_fish_sel(size))*(1-discard_survival) + temp_catch_mat(size)*(1-retain_fish_sel(size))*(1-discard_survival);
	   temp_imm(size) = temp_imm(size) * exp(-(exp(log_f+f_dev(year))*total_fish_sel(size)));
	   temp_mat(size) = temp_mat(size) * exp(-(exp(log_f+f_dev(year))*total_fish_sel(size)));
	   temp_imm(size) += temp_catch_imm(size)*(1-retain_fish_sel(size))*(discard_survival)	;
	   temp_mat(size) += temp_catch_mat(size)*(1-retain_fish_sel(size))*(discard_survival);
	   }	
	   for (int size=1;size<=size_n;size++) 
	   {
       temp_noskip(size) = temp_imm(size) * (1-prop_skip(year,size));
	   temp_skip(size) = temp_imm(size) * (prop_skip(year,size));
	   }
	  // growth
	   trans_imm = size_trans * temp_noskip;
	   // recruitment
       trans_imm(1) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(1);
       trans_imm(2) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(2);
	   trans_imm(3) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(3);
	  // maturity
	   for (int size=1;size<=size_n;size++) 
	   {
	  	temp_imm(size) = temp_skip(size) + (trans_imm(size) * (1-use_term_molt(year,size))); 	
	  	temp_mat(size) = (trans_imm(size) * use_term_molt(year,size) + temp_mat(size));
       }
	   // natural mortality		
	  for (int size=1;size<=size_n;size++) 
	   {
	   imm_n_size_pred(year+1,size) = temp_imm(size) * exp(-1*(0.41)*exp(nat_m(year,size)));
	   mat_n_size_pred(year+1,size) = temp_mat(size) * exp(-1*(0.41)*exp(nat_m_mat(year,size)));
       }
    }
   // cout<<"imm_n_size_pred"<<imm_n_size_pred<<endl;
   // cout<<"mat_n_size_pred"<<imm_n_size_pred<<endl;
   // cout<<"pred_retained_size_comp"<<pred_retained_size_comp<<endl;
   // cout<<"pred_discard_size_comp"<<pred_discard_size_comp<<endl;
   evaluate_the_objective_function();
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::evaluate_the_objective_function(void)
{
  // make total numbers by maturity state from obs and preds
  imm_numbers_pred.initialize();
  mat_numbers_pred.initialize();
  sum_imm_numbers_obs.initialize();
  sum_mat_numbers_obs.initialize(); 
  pred_retained_n.initialize();
  pred_discard_n.initialize();
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
   {
    imm_numbers_pred(year)    += selectivity(year,size)*imm_n_size_pred(year,size);
    mat_numbers_pred(year)    += selectivity_mat(year,size)*mat_n_size_pred(year,size);
	sum_imm_numbers_obs(year) += imm_n_size_obs(year,size);
	sum_mat_numbers_obs(year) += mat_n_size_obs(year,size);
   }
   for (int year=styr;year<endyr;year++)
   for (int size=1;size<=size_n;size++)
   { 
	pred_retained_n(year)     += pred_retained_size_comp(year,size);
	pred_discard_n(year)      += pred_discard_size_comp(year,size);
	   }
    // cout<<"imm_numbers_pred"<<imm_numbers_pred<<endl;
    // cout<<"imm_n_obs"<<imm_n_obs<<endl;
    // cout<<"mat_numbers_pred"<<mat_numbers_pred<<endl;
    // cout<<"mat_n_obs"<<mat_n_obs<<endl; 
    // cout<<"pred_retained_n"<<pred_retained_n<<endl;
    // cout<<"ret_cat_numbers"<<ret_cat_numbers<<endl;
    // cout<<"pred_discard_n"<<pred_discard_n<<endl;
    // cout<<"disc_cat_numbers"<<disc_cat_numbers<<endl;  
  // likelihoods
  imm_num_like = 0;
  for (int year=styr;year<=endyr;year++)
    imm_num_like += square( log(imm_numbers_pred(year)) - log(imm_n_obs(year))) / (2.0 * square(sigma_numbers_imm));
  mat_num_like = 0;
  for (int year=styr;year<=endyr;year++)
    mat_num_like += square( log(mat_numbers_pred(year)) - log(mat_n_obs(year))) / (2.0 * square(sigma_numbers_mat));
   // likelihoods
  ret_cat_like = 0;
  for (int year=cat_st;year<=d_cat_end;year++)
    ret_cat_like += square( log(pred_retained_n(year)) - log(ret_cat_numbers(year))) / (2.0 * square(sigma_numbers_ret));
  disc_cat_like = 0;
  for (int year=d_cat_st;year<=d_cat_end;year++)
    disc_cat_like += square( log(pred_discard_n(year)) - log(disc_cat_numbers(year))) / (2.0 * square(sigma_numbers_disc));
    // cout<<"imm_n_size_obs"<<imm_n_size_obs<<endl;
    // cout<<"imm_n_size_pred"<<imm_n_size_pred<<endl;
    // cout<<"mat_n_size_obs"<<mat_n_size_obs<<endl;
    // cout<<"mat_n_size_pred"<<mat_n_size_pred<<endl; 
    // cout<<"ret_cat_size"<<ret_cat_size<<endl;
    // cout<<"pred_retained_size_comp"<<pred_retained_size_comp<<endl;
    // cout<<"disc_cat_size"<<disc_cat_size<<endl;
    // cout<<"pred_discard_n"<<pred_discard_n<<endl;  
  // immature numbers at size data
   imm_like = 0;
   for (int year=styr;year<=endyr;year++)
    for (int size=1;size<=size_n;size++)
     if (imm_n_size_obs(year,size) >0 & imm_n_size_pred(year,size) >0)
      imm_like += imm_eff_samp*(imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)) * log( (selectivity(year,size)*imm_n_size_pred(year,size)/imm_numbers_pred(year)) / (imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)));
   imm_like = -1*imm_like;
  // mature numbers at size data
   mat_like = 0;
   for (int year=styr;year<=endyr;year++)
    for (int size=1;size<=size_n;size++)
     if (mat_n_size_obs(year,size) >0)
      mat_like += mat_eff_samp*(mat_n_size_obs(year,size)/sum_mat_numbers_obs(year)) * log( (selectivity_mat(year,size)*mat_n_size_pred(year,size)/mat_numbers_pred(year)) / (mat_n_size_obs(year,size)/sum_mat_numbers_obs(year)));
   mat_like = -1*mat_like;
  // retained catch at size data
   ret_comp_like = 0;
   for (int year=d_cat_st;year<d_cat_end;year++)
    for (int size=1;size<=size_n;size++)
     if (ret_cat_size(year,size) >0.001 & pred_retained_size_comp(year,size) >0.001)
      ret_comp_like += ret_eff_samp*(ret_cat_size(year,size)) * log( (pred_retained_size_comp(year,size)/pred_retained_n(year)) / (ret_cat_size(year,size)));
   ret_comp_like = -1*ret_comp_like;
  // discard catch at size data
   disc_comp_like = 0;
   for (int year=d_cat_st;year<d_cat_end;year++)
    for (int size=1;size<=size_n;size++)
	  if (disc_cat_size(year,size) >0.001 & pred_discard_size_comp(year,size) >0.001)
      disc_comp_like += ret_eff_samp*(disc_cat_size(year,size)) * log( (pred_discard_size_comp(year,size)/pred_discard_n(year)) / (disc_cat_size(year,size)));
   disc_comp_like = -1*disc_comp_like;
 //penalties on m 
  nat_m_mu_like =0;
  nat_m_mu_like += pow(((log_m_mu)-(log_mu_m_prior))/ (sqrt(2)*sqrt(sigma_m_mu(1))),2.0);
  nat_m_mat_mu_like =0;
  nat_m_mat_mu_like += pow(((log_m_mat_mu(1))-(log_mu_m_mat_prior(1)))/ (sqrt(2)*sqrt(sigma_m_mu(2))),2.0); 
  if(est_m_devs>0)
  {
  nat_m_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_like += pow((nat_m(year,1)-log_m_mu)/ (sqrt(2)*sqrt(sigma_m(1))),2.0);
  nat_m_mat_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_mat_like += pow((nat_m_mat(year,1)-log_m_mat_mu(1))/ (sqrt(2)*sqrt(sigma_m(2))),2.0);
  nat_m_lg_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_lg_like += pow((nat_m_mat(year,18)-log_m_mat_mu(2))/ (sqrt(2)*sqrt(sigma_m(3))),2.0);
  }
  if(est_q_devs>0)
  {
  q_like =0;
  for (int year=styr;year<=endyr;year++)
   q_like += pow((selectivity(year,4)-surv_sel(4))/ (sqrt(2)*sqrt(sigma_q(1))),2.0);
  q_mat_like =0;
  for (int year=styr;year<=endyr;year++)
   q_mat_like += pow((selectivity_mat(year,4)-surv_sel(4))/ (sqrt(2)*sqrt(sigma_q(2))),2.0);
  }
  smooth_q_like = 0;
  smooth_q_like = smooth_q_weight* (norm2(first_difference(first_difference(q_dev))) +norm2(first_difference(first_difference(q_mat_dev)))) ;
  smooth_m_like = 0;
  smooth_m_like = smooth_m_weight(1)* norm2(first_difference(first_difference(nat_m_dev))) + smooth_m_weight(2)*norm2(first_difference(first_difference(nat_m_mat_dev))) ;
  smooth_f_like = 0;
  smooth_f_like = smooth_f_weight* (norm2(first_difference(first_difference(f_dev))));
  smooth_surv_like = 0;
  smooth_surv_like = smooth_surv_weight* (norm2(first_difference(first_difference(surv_sel))));
  f = imm_num_like + mat_num_like + ret_cat_like + disc_cat_like + imm_like + mat_like + ret_comp_like + disc_comp_like + 
  nat_m_like + nat_m_mat_like + nat_m_lg_like + nat_m_mu_like + nat_m_mat_mu_like + smooth_q_like + smooth_m_like + q_like + q_mat_like + smooth_f_like +
  surv_sel_prior + smooth_surv_like;
  cout<<imm_num_like<< " " << mat_num_like << " " << imm_like << " " << mat_like << " " <<endl;
  cout<<ret_cat_like<< " " << ret_cat_like << " " << ret_comp_like << " " << disc_comp_like << " " <<endl;
  cout<<nat_m_like<< " " << disc_cat_like << " " << nat_m_lg_like << " " << nat_m_mu_like << " " << nat_m_mat_mu_like << " " <<endl;
  cout<<smooth_q_like<< " " << smooth_m_like << " " << q_like << " " << q_mat_like << " " <<smooth_f_like<<" " << surv_sel_prior<<" "<<endl;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report<<"$likelihoods"<<endl;
  report<<f<<" "<<imm_num_like<<" "<<mat_num_like<<" "<<imm_like<<" "<<mat_like<<" "<<nat_m_like<<" "<<nat_m_mat_like<<" "<<nat_m_lg_like<< " " << nat_m_mu_like<<" "<<nat_m_mat_mu_like<<" "<<smooth_q_like<<" "<<smooth_m_like<<" "<<q_like<<" "<<q_mat_like<<endl;
  report <<"$natural mortality" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(nat_m(i))<<endl;
  }
  report <<"$size_trans" << endl;
  for(int i=1; i<=size_n; i++)
  {
    report << size_trans(i)<<endl;
  }
  report <<"$recruits" << endl;
  for(int i=styr; i<=endyr; i++)
   report << mfexp(log_avg_rec + rec_devs(i))<<endl;
  report <<"$immature numbers at size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (elem_prod(selectivity(i),imm_n_size_pred(i)))/imm_numbers_pred(i)<<endl;
  }
  report <<"$mature numbers at size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (elem_prod(selectivity_mat(i),mat_n_size_pred(i)))/mat_numbers_pred(i)<<endl;
  }
  report<<"$styr"<<endl;
  report<<styr<<endl;
  report<<"$endyr"<<endl;
  report<<endyr<<endl;
  report<<"$imm_numbers_pred"<<endl;
  report<<imm_numbers_pred<<endl;
  report<<"$mat_numbers_pred"<<endl;
  report<<mat_numbers_pred<<endl;
  report<<"$imm_n_obs"<<endl;
  report<<imm_n_obs<<endl;
  report<<"$mat_n_obs"<<endl;
  report<<mat_n_obs<<endl;
 // report<<"$survey_sel_input"<<endl;
  //report<<surv_sel<<endl;
  report<<"$ret_fish_sel"<<endl;
  report<<retain_fish_sel<<endl;
  report<<"$total_fish_sel"<<endl;
  report<<total_fish_sel<<endl;  
  report<<"$surv_sel"<<endl;
  report<<surv_sel<<endl;  
  report <<"$est_fishing_mort" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(log_f + f_dev(i))<<endl;
  }  
  report <<"$mature natural mortality" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(nat_m_mat(i))<<endl;
  }
  report <<"$survey selectivity" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << selectivity(i)<<endl;
  }
  report <<"$mature survey selectivity" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << selectivity_mat(i)<<endl;
  }
  report <<"$pred_imm_pop_num" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (imm_n_size_pred(i))<<endl;
  }
  report <<"$pred_mat_pop_num" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (mat_n_size_pred(i))<<endl;
  }
  report <<"$pred_retained_n" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (pred_retained_n(i))<<endl;
  }
  report <<"$pred_discard_n" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (pred_discard_n(i))<<endl;
  }
  report <<"$ret_cat_numbers" << endl;
  for(int i=styr; i<=d_cat_end; i++)
  {
    report << (ret_cat_numbers(i))<<endl;
  }
  report <<"$disc_cat_numbers" << endl;
  for(int i=d_cat_st; i<=d_cat_end; i++)
  {
    report << (disc_cat_numbers(i))<<endl;
  }
  report <<"$pred_retained_size_comp" << endl;
  for(int i=styr; i<=d_cat_end; i++)
  {
    report << ((pred_retained_size_comp(i)/pred_retained_n(i)))<<endl;
  }
  report <<"$pred_discard_size_comp" << endl;
  for(int i=styr; i<=d_cat_end; i++)
  {
    report << ((pred_discard_size_comp(i)/pred_discard_n(i)))<<endl;
  }
  report <<"$obs_retained_size_comp" << endl;
  for(int i=d_cat_st; i<=d_cat_end; i++)
  {
    report << (ret_cat_size(i))<<endl;
  }
  report <<"$obs_discard_size_comp" << endl;
  for(int i=d_cat_st; i<=d_cat_end; i++)
  {
    report << (disc_cat_size(i))<<endl;
  }
  report <<"$temp_prop_rec" << endl;
  report << temp_prop_rec << endl;
  report <<"$prob_term_molt" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (use_term_molt(i))<<endl;
  }
    save_gradients(gradients);
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
