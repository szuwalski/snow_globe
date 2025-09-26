DATA_SECTION
 //Bering Sea snow crab model
  
  init_int styr  																							
  init_int endyr 
  init_int dat_yr
  init_ivector years(1,dat_yr)
  init_int size_n
  init_vector sizes(1,size_n)
  init_vector imm_n_obs(styr,endyr)
  init_matrix imm_n_size_obs(styr,endyr,1,size_n)
  init_vector mat_n_obs(styr,endyr)
  init_matrix mat_n_size_obs(styr,endyr,1,size_n)
  init_matrix prop_term_molt(styr,endyr,1,size_n)
  init_matrix size_trans(1,size_n,1,size_n) 
  //init_vector survey_sel(1,size_n)
  init_matrix skip_molt(styr,endyr,1,size_n) 
  init_number sigma_numbers_imm
  init_number sigma_numbers_mat

  init_number mat_eff_samp
  init_number imm_eff_samp

  init_vector log_mu_m_prior(1,3)
  init_number est_m_devs
  init_number est_q_devs
  init_number est_m_mat_devs
  init_number est_q_mat_devs
  init_number est_sigma_m
  init_vector sigma_m_mu(1,3)
  init_number smooth_q_weight
  init_vector smooth_m_weight(1,2)
  init_number est_log_m_mu
  init_number est_sigma_q
  init_number est_m_lg_devs
  init_number smooth_f_weight
 // init_number surv_sel_cv
  init_number smooth_surv_weight
  init_number est_sel
   // !!cout<<"imm_n_obs"<<imm_n_obs<<endl; 
  // !!cout<<"mat_n_obs"<<mat_n_obs<<endl; 
  // !!cout<<"est_m_devs"<<est_m_devs<<endl;
    // !!cout<<"sigma_m_mu"<<sigma_m_mu<<endl;
	  // !!cout<<"smooth_m_weight"<<smooth_m_weight<<endl;
	     !!cout<<"smooth_surv_weight"<<smooth_surv_weight<<endl;

 //==read in catch data
 !! ad_comm::change_datafile_name("catch_dat.DAT");
  init_int ret_cat_yr_n
  init_int cat_st
  init_int cat_end
  init_vector ret_cat_yrs(1,ret_cat_yr_n)
  init_vector ret_cat_numbers(cat_st,cat_end)
   
  init_int disc_cat_yr_n
  init_int d_cat_st
  init_int d_cat_end 
  init_vector disc_cat_yrs(1,disc_cat_yr_n)
  init_vector disc_cat_numbers(d_cat_st,d_cat_end) 
  
  init_number sigma_numbers_ret
  init_number sigma_numbers_disc
  init_number discard_survival
  init_number ret_eff_samp
  init_number disc_eff_samp
  
  //uses discards because observers only for size comps
  init_matrix ret_cat_size(d_cat_st,d_cat_end,1,size_n)
  init_matrix disc_cat_size(d_cat_st,d_cat_end,1,size_n)
  
      !!cout<<"sigma_numbers_ret"<<sigma_numbers_ret<<endl;
	  !!cout<<"disc_eff_samp"<<disc_eff_samp<<endl;
  
PARAMETER_SECTION
  init_bounded_vector log_n_imm(1,size_n,1.01,30,1)
  init_bounded_vector log_n_mat(1,size_n,1.01,30,1)
  init_bounded_dev_vector nat_m_dev(styr,endyr,-4,4,est_m_devs)
  init_bounded_dev_vector nat_m_mat_dev(styr,endyr,-4,4,est_m_mat_devs)
  init_bounded_dev_vector nat_m_lg_dev(styr,endyr,-4,4,est_m_mat_devs)
  init_bounded_vector q_dev(styr,endyr,-0.2,0.2,est_q_devs)
  init_bounded_vector q_mat_dev(styr,endyr,-0.2,0.2,est_q_mat_devs)
  init_bounded_vector q_lg_dev(styr,endyr,-0.2,0.2,est_q_mat_devs)
  init_bounded_number log_avg_rec(1,40)
  init_bounded_dev_vector rec_devs(styr,endyr,-10,10,1)
  init_bounded_vector sigma_m(1,3,0.01,4,est_sigma_m)
  init_bounded_vector log_m_mu(1,3,-5,3,est_log_m_mu)
  init_bounded_vector prop_rec(1,2,0.00001,200)
  init_bounded_vector sigma_q(1,2,0.01,4,est_sigma_q)
  
  init_bounded_number log_f(-5,5)
  init_bounded_dev_vector f_dev(styr,endyr,-5,5)
  init_bounded_number fish_ret_sel_50(25,150,2)
  init_bounded_number fish_ret_sel_slope(0.0001,20,-2)
  init_bounded_number fish_tot_sel_50(25,150,2)
  init_bounded_number fish_tot_sel_slope(0.0001,20,2) 
  //init_bounded_vector surv_sel(1,size_n,0.000001,1)

  init_bounded_number surv_omega(0,1,est_sel)
  init_bounded_number surv_alpha1(0.0001,20,est_sel)
  init_bounded_number surv_beta1(25,150,est_sel)
  init_bounded_number surv_alpha2(0.0001,20,est_sel)
  init_bounded_number surv_beta2(25,150,est_sel)

  matrix imm_n_size_pred(styr,endyr,1,size_n)
  matrix mat_n_size_pred(styr,endyr,1,size_n)
  matrix nat_m(styr,endyr,1,size_n)
  matrix nat_m_mat(styr,endyr,1,size_n)
  matrix selectivity(styr,endyr,1,size_n)
  matrix selectivity_mat(styr,endyr,1,size_n)
  
  vector total_fish_sel(1,size_n)
  vector retain_fish_sel(1,size_n)
  vector pred_retained_n(styr,endyr)
  vector pred_discard_n(styr,endyr)
  matrix pred_retained_size_comp(styr,endyr-1,1,size_n)
  matrix pred_discard_size_comp(styr,endyr-1,1,size_n)
  
  vector surv_sel(1,size_n)
  vector temp_imm(1,size_n)
  vector temp_mat(1,size_n)
  vector trans_imm(1,size_n)
  vector skip_imm(1,size_n)
  vector temp_catch_imm(1,size_n)
  vector temp_catch_mat(1,size_n)

  vector sum_imm_numbers_obs(styr,endyr)
  vector sum_mat_numbers_obs(styr,endyr)
  vector imm_numbers_pred(styr,endyr)
  vector mat_numbers_pred(styr,endyr)
  vector sum_ret_numbers_obs(styr,endyr)
  vector sum_disc_numbers_obs(styr,endyr) 
  
  number imm_num_like
  number mat_num_like
  number ret_cat_like
  number disc_cat_like
  number imm_like
  number mat_like
  number ret_comp_like
  number disc_comp_like
  
  matrix use_term_molt(styr,endyr,1,size_n)
  
  number nat_m_like
  number nat_m_mat_like
  number nat_m_mu_like
  number nat_m_lg_like
  number nat_m_mat_mu_like
  number nat_m_lg_mu_like
  number smooth_q_like
  number smooth_m_like
  number q_like
  number q_mat_like
  number smooth_f_like
  number surv_sel_prior
  number smooth_surv_like
  
  vector temp_prop_rec(1,3)
  number tot_prop_rec
  
  sdreport_vector total_population_n(styr,endyr)
  sdreport_vector fished_population_n(styr,endyr) 
  sdreport_vector recruits(styr,endyr)  
  
  objective_function_value f
  
 
 
//==============================================================================
PROCEDURE_SECTION
 dvariable fpen=0.0;
// initial year numbers at size and selectivity
 for(int size=1;size<=size_n;size++)
  {
   imm_n_size_pred(styr,size) = exp(log_n_imm(size));
   mat_n_size_pred(styr,size) = exp(log_n_mat(size));
   total_fish_sel(size) = 1 / (1+exp(-fish_tot_sel_slope*(sizes(size)-fish_tot_sel_50))) ; 
   retain_fish_sel(size) = 1 / (1+exp(-fish_ret_sel_slope*(sizes(size)-fish_ret_sel_50)));
   surv_sel(size) = (surv_omega / (1 + exp(-surv_alpha1*(sizes(size)-surv_beta1))) ) + ( (1-surv_omega)/(1+exp(-surv_alpha2*(sizes(size)-surv_beta2)))  );
   }

 for(int year=styr;year<=endyr;year++)
  for(int size=1;size<=size_n;size++)
  {
  nat_m(year,size) = log_m_mu(1) + nat_m_dev(year);
  nat_m_mat(year,size) = log_m_mu(2) + nat_m_dev(year);

  selectivity(year,size) = surv_sel(size) + q_dev(year);
  selectivity_mat(year,size) = surv_sel(size) + q_dev(year);
  
  //==option for estimating separate devs by maturity state
  if(est_q_mat_devs>0)
   selectivity_mat(year,size) = surv_sel(size) + q_mat_dev(year);  
  if(est_m_mat_devs>0)
   nat_m_mat(year,size) = log_m_mu(2) + nat_m_mat_dev(year);
  if(est_m_lg_devs>0 & size > 13)
  {
   nat_m(year,size) = log_m_mu(3) + nat_m_lg_dev(year);
   nat_m_mat(year,size) = log_m_mu(3) + nat_m_lg_dev(year);
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
		temp_imm(size) = imm_n_size_pred(year,size) * exp(-1*(0.67)*exp( nat_m(year,size)));
		temp_mat(size) = mat_n_size_pred(year,size) * exp(-1*(0.67)*exp( nat_m_mat(year,size)));
	   }	

  skip_imm.initialize();
  for (int size=1;size<=size_n;size++) 
      {
      skip_imm(size) = skip_molt(year,size)*temp_imm(size);
	  temp_imm(size) = (1-skip_molt(year,size))*temp_imm(size);
	  }

	  // growth
	   trans_imm = size_trans * temp_imm;
	   
	   // recruitment

       trans_imm(1) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(1);
       trans_imm(2) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(2);
	   trans_imm(3) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(3);
       recruits(year) = exp(log_avg_rec + rec_devs(year));
	  // maturity
	   for (int size=1;size<=size_n;size++) 
	   {
	  	temp_imm(size) = (trans_imm(size) * (1-use_term_molt(year,size))) + skip_imm(size); 	
	  	temp_mat(size) = (trans_imm(size) * use_term_molt(year,size) + temp_mat(size));
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

	   // natural mortality		
	  for (int size=1;size<=size_n;size++) 
	   {
	   imm_n_size_pred(year+1,size) = temp_imm(size) * exp(-1*(0.33)*exp(nat_m(year,size)));
	   mat_n_size_pred(year+1,size) = temp_mat(size) * exp(-1*(0.33)*exp(nat_m_mat(year,size)));
       }
	   

    }

   // cout<<"imm_n_size_pred"<<imm_n_size_pred<<endl;
   // cout<<"mat_n_size_pred"<<imm_n_size_pred<<endl;
   // cout<<"pred_retained_size_comp"<<pred_retained_size_comp<<endl;
   // cout<<"pred_discard_size_comp"<<pred_discard_size_comp<<endl;
   evaluate_the_objective_function();
   
//==============================================================================
FUNCTION evaluate_the_objective_function

  // make total numbers by maturity state from obs and preds
  imm_numbers_pred.initialize();
  mat_numbers_pred.initialize();
  sum_imm_numbers_obs.initialize();
  sum_mat_numbers_obs.initialize(); 
  pred_retained_n.initialize();
  pred_discard_n.initialize();
  total_population_n.initialize();
  fished_population_n.initialize();
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
   {
    imm_numbers_pred(year)    += selectivity(year,size)*imm_n_size_pred(year,size);
    mat_numbers_pred(year)    += selectivity_mat(year,size)*mat_n_size_pred(year,size);
	sum_imm_numbers_obs(year) += imm_n_size_obs(year,size);
	sum_mat_numbers_obs(year) += mat_n_size_obs(year,size);
		total_population_n(year)  += imm_n_size_pred(year,size)+mat_n_size_pred(year,size);
	fished_population_n(year)  += retain_fish_sel(size)*imm_n_size_pred(year,size)+retain_fish_sel(size)*mat_n_size_pred(year,size);

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
  nat_m_mu_like += pow(((log_m_mu(1))-(log_mu_m_prior(1)))/ (sqrt(2)*sqrt(sigma_m_mu(1))),2.0);
  nat_m_mat_mu_like =0;
  nat_m_mat_mu_like += pow(((log_m_mu(2))-(log_mu_m_prior(2)))/ (sqrt(2)*sqrt(sigma_m_mu(2))),2.0); 
   
  if(est_m_devs>0 & current_phase()>=est_m_devs)
  {
  nat_m_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_like += pow((nat_m(year,1)-log_m_mu(1))/ (sqrt(2)*sqrt(sigma_m(1))),2.0);
   
  nat_m_mat_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_mat_like += pow((nat_m_mat(year,1)-log_m_mu(2))/ (sqrt(2)*sqrt(sigma_m(2))),2.0);

  nat_m_lg_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_lg_like += pow((nat_m_mat(year,20)-log_m_mu(3))/ (sqrt(2)*sqrt(sigma_m(3))),2.0);
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
  
// ========================y==================================================   
REPORT_SECTION
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
  
  report <<"$skip_molt" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (skip_molt(i))<<endl;
  }
  
    report <<"$obs_imm_n_size" << endl;
  for(int i=styr; i<=endyr; i++)
  {

    report << imm_n_size_obs(i)<<endl;
  }
  
  report <<"$obs_mat_n_size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mat_n_size_obs(i)<<endl;
  }

  
   report <<"$sizes" << endl;
  report << sizes << endl;	
  
  report <<"$imm_cv" << endl;
  report << sigma_numbers_imm << endl;	
  
  report <<"$mat_cv" << endl;
  report << sigma_numbers_mat << endl;	
 
  report <<"$ret_cat_yrs" << endl;
  report << ret_cat_yrs << endl;	 
  	 
     report <<"$disc_cat_yrs" << endl;
  report << disc_cat_yrs << endl;	 
    save_gradients(gradients);
  
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 10000
  convergence_criteria 1e-3

