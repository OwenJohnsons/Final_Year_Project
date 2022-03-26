 #!/bin/csh
  
 if (! -e Geneva) mkdir Geneva
 cp -p ../Data/Geneva/geneva.iso                                                               Geneva/geneva.iso                                                               
 if (! -e DAM97) mkdir DAM97
 cp -p ../Data/DAM97/dam97.dat                                                                 DAM97/dam97.dat                                                                 
 if (! -e Lyon) mkdir Lyon
 cp -p ../Data/Lyon/baraffe98alpha1.0.iso                                                      Lyon/baraffe98alpha1.0.iso                                                      
 if (! -e Lyon) mkdir Lyon
 cp -p ../Data/Lyon/baraffe98alpha1.9.iso                                                      Lyon/baraffe98alpha1.9.iso                                                      
 cp -p ../Data/palla_stahler.dat                                                               palla_stahler.dat                                                               
 if (! -e Lyon) mkdir Lyon
 cp -p ../Data/Lyon/dusty00.iso                                                                Lyon/dusty00.iso                                                                
 if (! -e Siess) mkdir Siess
 cp -p ../Data/Siess/siessz02.dat                                                              Siess/siessz02.dat                                                              
 if (! -e Siess) mkdir Siess
 cp -p ../Data/Siess/siessz01.dat                                                              Siess/siessz01.dat                                                              
 if (! -e Padova) mkdir Padova
 cp -p ../Data/Padova/isoc_z019_new.iso                                                        Padova/isoc_z019_new.iso                                                        
 if (! -e Lyon) mkdir Lyon
 cp -p ../Data/Lyon/cond03.iso                                                                 Lyon/cond03.iso                                                                 
 cp -p ../Data/dartmouth_z0.018.trk                                                            dartmouth_z0.018.trk                                                            
 cp -p ../Data/dartmouth_z0.026.trk                                                            dartmouth_z0.026.trk                                                            
 cp -p ../Data/brott.iso                                                                       brott.iso                                                                       
 cp -p ../Data/pisa_z0.013.dat                                                                 pisa_z0.013.dat                                                                 
 cp -p ../Data/tycho_bessell_Mp0.0.bc                                                          tycho_bessell_Mp0.0.bc                                                          
 cp -p ../Data/tycho.rv                                                                        tycho.rv                                                                        
 cp -p ../Data/kenyon_hartmann_Mp0.0.bc                                                        kenyon_hartmann_Mp0.0.bc                                                        
 cp -p ../Data/with_isochrones.rv                                                              with_isochrones.rv                                                              
 cp -p ../Data/btsettl_Mp0.0_bessell_energy.bc                                                 btsettl_Mp0.0_bessell_energy.bc                                                 
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_bessell_energy.ext                                     Ext_Tables/btsettl_Mp0.0_bessell_energy.ext                                     
 cp -p ../Data/btsettl_Mp0.0_bessell_photon.bc                                                 btsettl_Mp0.0_bessell_photon.bc                                                 
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     
 cp -p ../Data/dartmouth_btsettl_Mp0.0_johnson_tuned_photon.bc                                 dartmouth_btsettl_Mp0.0_johnson_tuned_photon.bc                                 
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     
 cp -p ../Data/newbaraffe_19_btsettl_Mp0.0_johnson_tuned_photon.bc                             newbaraffe_19_btsettl_Mp0.0_johnson_tuned_photon.bc                             
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     
 cp -p ../Data/pisa_btsettl_Mp0.0_johnson_tuned_photon.bc                                      pisa_btsettl_Mp0.0_johnson_tuned_photon.bc                                      
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     
 cp -p ../Data/btsettl_Mp0.0_wfc.bc                                                            btsettl_Mp0.0_wfc.bc                                                            
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_wfc.ext                                                Ext_Tables/btsettl_Mp0.0_wfc.ext                                                
 cp -p ../Data/dam97_btsettl_Mp0.0_wfc_tuned_photon.bc                                         dam97_btsettl_Mp0.0_wfc_tuned_photon.bc                                         
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_wfc.ext                                                Ext_Tables/btsettl_Mp0.0_wfc.ext                                                
 cp -p ../Data/dartmouth_btsettl_Mp0.0_wfc_tuned_photon.bc                                     dartmouth_btsettl_Mp0.0_wfc_tuned_photon.bc                                     
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_wfc.ext                                                Ext_Tables/btsettl_Mp0.0_wfc.ext                                                
 cp -p ../Data/newbaraffe_19_btsettl_Mp0.0_wfc_tuned_photon.bc                                 newbaraffe_19_btsettl_Mp0.0_wfc_tuned_photon.bc                                 
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_wfc.ext                                                Ext_Tables/btsettl_Mp0.0_wfc.ext                                                
 cp -p ../Data/pisa_btsettl_Mp0.0_wfc_tuned_photon.bc                                          pisa_btsettl_Mp0.0_wfc_tuned_photon.bc                                          
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_wfc.ext                                                Ext_Tables/btsettl_Mp0.0_wfc.ext                                                
 cp -p ../Data/btsettl_Mp0.0_iphas_uvex.bc                                                     btsettl_Mp0.0_iphas_uvex.bc                                                     
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         
 cp -p ../Data/dam97_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                                  dam97_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                                  
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         
 cp -p ../Data/dartmouth_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                              dartmouth_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                              
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         
 cp -p ../Data/newbaraffe_19_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                          newbaraffe_19_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                          
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         
 cp -p ../Data/pisa_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                                   pisa_btsettl_Mp0.0_iphas_uvex_tuned_photon.bc                                   
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         Ext_Tables/btsettl_Mp0.0_iphas_uvex.ext                                         
 cp -p ../Data/btsettl_Mp0.0_sdss.bc                                                           btsettl_Mp0.0_sdss.bc                                                           
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_sdss.ext                                               Ext_Tables/btsettl_Mp0.0_sdss.ext                                               
 cp -p ../Data/baraffe_btsettl_Mp0.0_sdss_tuned.bc                                             baraffe_btsettl_Mp0.0_sdss_tuned.bc                                             
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_sdss.ext                                               Ext_Tables/btsettl_Mp0.0_sdss.ext                                               
 cp -p ../Data/dartmouth_btsettl_Mp0.0_sdss_tuned.bc                                           dartmouth_btsettl_Mp0.0_sdss_tuned.bc                                           
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_sdss.ext                                               Ext_Tables/btsettl_Mp0.0_sdss.ext                                               
 cp -p ../Data/pisa_btsettl_Mp0.0_sdss_tuned.bc                                                pisa_btsettl_Mp0.0_sdss_tuned.bc                                                
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_sdss.ext                                               Ext_Tables/btsettl_Mp0.0_sdss.ext                                               
 cp -p ../Data/pisa_btsettl_Mp0.0_johnson_tuned_photon.bc                                      pisa_btsettl_Mp0.0_johnson_tuned_photon.bc                                      
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     Ext_Tables/btsettl_Mp0.0_bessell_photon.ext                                     
 cp -p ../Data/btsettl_Mp0.0_cfht.bc                                                           btsettl_Mp0.0_cfht.bc                                                           
 cp -p ../Data/                                                                                                                                                                
 cp -p ../Data/btsettl_Mp0.0_ukidss.bc                                                         btsettl_Mp0.0_ukidss.bc                                                         
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_ukidss.ext                                             Ext_Tables/btsettl_Mp0.0_ukidss.ext                                             
 cp -p ../Data/btsettl_Mp0.0_2mass.bc                                                          btsettl_Mp0.0_2mass.bc                                                          
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_2mass.ext                                              Ext_Tables/btsettl_Mp0.0_2mass.ext                                              
 cp -p ../Data/dartmouth_btsettl_Mp0.0_2mass_tuned.bc                                          dartmouth_btsettl_Mp0.0_2mass_tuned.bc                                          
 cp -p ../Data/rieke_lebofsky.rv                                                               rieke_lebofsky.rv                                                               
 cp -p ../Data/newbaraffe_19_btsettl_Mp0.0_2mass_tuned.bc                                      newbaraffe_19_btsettl_Mp0.0_2mass_tuned.bc                                      
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_2mass.ext                                              Ext_Tables/btsettl_Mp0.0_2mass.ext                                              
 cp -p ../Data/pisa_btsettl_Mp0.0_2mass_tuned.bc                                               pisa_btsettl_Mp0.0_2mass_tuned.bc                                               
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/btsettl_Mp0.0_2mass.ext                                              Ext_Tables/btsettl_Mp0.0_2mass.ext                                              
 cp -p ../Data/kurucz_Mp0.0_bessell_energy.bc                                                  kurucz_Mp0.0_bessell_energy.bc                                                  
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/kurucz_Mp0.0_bessell_energy.ext                                      Ext_Tables/kurucz_Mp0.0_bessell_energy.ext                                      
 cp -p ../Data/kurucz_Mp0.0_bessell_photon.bc                                                  kurucz_Mp0.0_bessell_photon.bc                                                  
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/kurucz_Mp0.0_bessell_photon.ext                                      Ext_Tables/kurucz_Mp0.0_bessell_photon.ext                                      
 cp -p ../Data/kurucz_Mp0.0_wfc.bc                                                             kurucz_Mp0.0_wfc.bc                                                             
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/kurucz_Mp0.0_wfc.ext                                                 Ext_Tables/kurucz_Mp0.0_wfc.ext                                                 
 cp -p ../Data/kurucz_Mp0.0_sdss.bc                                                            kurucz_Mp0.0_sdss.bc                                                            
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/kurucz_Mp0.0_sdss.ext                                                Ext_Tables/kurucz_Mp0.0_sdss.ext                                                
 cp -p ../Data/kurucz_Mp0.0_ukidss.bc                                                          kurucz_Mp0.0_ukidss.bc                                                          
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/kurucz_Mp0.0_ukidss.ext                                              Ext_Tables/kurucz_Mp0.0_ukidss.ext                                              
 cp -p ../Data/kurucz_Mp0.0_2mass.bc                                                           kurucz_Mp0.0_2mass.bc                                                           
 if (! -e Ext_Tables) mkdir Ext_Tables
 cp -p ../Data/Ext_Tables/kurucz_Mp0.0_2mass.ext                                               Ext_Tables/kurucz_Mp0.0_2mass.ext                                               
