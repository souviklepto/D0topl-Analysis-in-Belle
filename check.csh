# /bin/csh -f 
#setenv USE_GRAND_REPROCESS_DATA 1
#source /sw/belle/local/etc/cshrc_general
#setenv BASF_USER_INIT geant_init
#setenv BASF_USER_IF basfsh.so
 
#setenv EVTGEN_PDL ${BELLE_TOP_DIR}/share/data-files/evtgenutil/evt.pdl
#setenv BELLE_USE_EVTGEN
 
   
 
#basf << EOF  >& bcs_sumchisq.log 
basf << EOF  #>& check.log 
module register user_ana_check
path create main
path add_module main fix_mdst user_ana_check
module put_parameter fix_mdst Endcap_MX_layer\11
module put_parameter fix_mdst Correct_ecl_pv\1
 
initialize 
 
#histogram define bcs_sumchisq.hbk
 histogram define check.hbk

#process_event /gpfs/home/belle2/souvik/PUBLIC/mcproduzh/gsim/mdst/evtgen_exp_09_dppippiz-380.mdst 

#process_event /gpfs/group/belle/bdata_b/dstprod/dat/e000009/HadronB/0416/on_resonance/00/e000009r000017-b20020416_1604.mdst

#process_event  /home/belle2/souvik/PUBLIC/mcproduzh/gsim/mdst/mdst/evtgen_exp_55_d0prel-0.mdst
#process_event  /group/belle/users/souvik/new_mode_mdst/new/new/kk_mode/evtgen_exp_37b_d0kk_1-0.mdst
process_event  /group/belle/users/souvik/new_mode_mdst/new/new/kk_mode/evtgen_exp_25_d0kk_1-17.mdst
#process_event  /group/belle/users/souvik/new_mode_mdst/new/new/muon_mode/evtgen_exp_37b_d0prmu_1-1.mdst

EOF
 
#exit 0
