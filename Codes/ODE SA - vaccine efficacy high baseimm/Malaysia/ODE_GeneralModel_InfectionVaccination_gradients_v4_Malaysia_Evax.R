COVID_vocBoost_hybrid_v4 <- function(t, x, parms) {

  #Uninfected
  Sch=x[1]; Scl=x[2]; Sah=x[3]; Sal=x[4]; Seh=x[5]; Sel=x[6] #S
  Ech=x[7]; Ecl=x[8]; Eah=x[9]; Eal=x[10]; Eeh=x[11]; Eel=x[12] #E
  Ach=x[13]; Acl=x[14]; Aah=x[15]; Aal=x[16]; Aeh=x[17]; Ael=x[18] #A
  Ich=x[19]; Icl=x[20]; Iah=x[21]; Ial=x[22]; Ieh=x[23]; Iel=x[24] #I
  Rch=x[25]; Rcl=x[26]; Rah=x[27]; Ral=x[28]; Reh=x[29]; Rel=x[30] #R
  Dch=x[31]; Dcl=x[32]; Dah=x[33]; Dal=x[34]; Deh=x[35]; Del=x[36] #D
  
  #Infected once, unvaccinated
  Sch2=x[37]; Scl2=x[38]; Sah2=x[39]; Sal2=x[40]; Seh2=x[41]; Sel2=x[42] #S
  Ech2=x[43]; Ecl2=x[44]; Eah2=x[45]; Eal2=x[46]; Eeh2=x[47]; Eel2=x[48] #E
  Ach2=x[49]; Acl2=x[50]; Aah2=x[51]; Aal2=x[52]; Aeh2=x[53]; Ael2=x[54] #A
  Ich2=x[55]; Icl2=x[56]; Iah2=x[57]; Ial2=x[58]; Ieh2=x[59]; Iel2=x[60] #I
  Rch2=x[61]; Rcl2=x[62]; Rah2=x[63]; Ral2=x[64]; Reh2=x[65]; Rel2=x[66] #R
  Dch2=x[67]; Dcl2=x[68]; Dah2=x[69]; Dal2=x[70]; Deh2=x[71]; Del2=x[72] #D  
  
  #Vaccinated with primary series, uninfected
  SVch=x[73]; SVcl=x[74]; SVah=x[75]; SVal=x[76]; SVeh=x[77]; SVel=x[78] #SV
  EVch=x[79]; EVcl=x[80]; EVah=x[81]; EVal=x[82]; EVeh=x[83]; EVel=x[84] #EV
  AVch=x[85]; AVcl=x[86]; AVah=x[87]; AVal=x[88]; AVeh=x[89]; AVel=x[90] #AV
  IVch=x[91]; IVcl=x[92]; IVah=x[93]; IVal=x[94]; IVeh=x[95]; IVel=x[96] #IV
  RVch=x[97]; RVcl=x[98]; RVah=x[99]; RVal=x[100]; RVeh=x[101]; RVel=x[102] #RV
  DVch=x[103]; DVcl=x[104]; DVah=x[105]; DVal=x[106]; DVeh=x[107]; DVel=x[108] #DV
  
  #Vaccinated with primary series and infected
  SVch2=x[109]; SVcl2=x[110]; SVah2=x[111]; SVal2=x[112]; SVeh2=x[113]; SVel2=x[114] #SV2
  EVch2=x[115]; EVcl2=x[116]; EVah2=x[117]; EVal2=x[118]; EVeh2=x[119]; EVel2=x[120] #EV2
  AVch2=x[121]; AVcl2=x[122]; AVah2=x[123]; AVal2=x[124]; AVeh2=x[125]; AVel2=x[126] #AV2
  IVch2=x[127]; IVcl2=x[128]; IVah2=x[129]; IVal2=x[130]; IVeh2=x[131]; IVel2=x[132] #IV2
  RVch2=x[133]; RVcl2=x[134]; RVah2=x[135]; RVal2=x[136]; RVeh2=x[137]; RVel2=x[138] #RV2
  DVch2=x[139]; DVcl2=x[140]; DVah2=x[141]; DVal2=x[142]; DVeh2=x[143]; DVel2=x[144] #DV2
  
  #Vaccinated with booster and uninfected
  SBch=x[145]; SBcl=x[146]; SBah=x[147]; SBal=x[148]; SBeh=x[149]; SBel=x[150] #SV
  EBch=x[151]; EBcl=x[152]; EBah=x[153]; EBal=x[154]; EBeh=x[155]; EBel=x[156] #EV
  ABch=x[157]; ABcl=x[158]; ABah=x[159]; ABal=x[160]; ABeh=x[161]; ABel=x[162] #AV
  IBch=x[163]; IBcl=x[164]; IBah=x[165]; IBal=x[166]; IBeh=x[167]; IBel=x[168] #IV
  RBch=x[169]; RBcl=x[170]; RBah=x[171]; RBal=x[172]; RBeh=x[173]; RBel=x[174] #RV
  DBch=x[175]; DBcl=x[176]; DBah=x[177]; DBal=x[178]; DBeh=x[179]; DBel=x[180] #DV
  
  SBch2=x[181]; SBcl2=x[182]; SBah2=x[183]; SBal2=x[184]; SBeh2=x[185]; SBel2=x[186] #SB2
  EBch2=x[187]; EBcl2=x[188]; EBah2=x[189]; EBal2=x[190]; EBeh2=x[191]; EBel2=x[192] #EB2
  ABch2=x[193]; ABcl2=x[194]; ABah2=x[195]; ABal2=x[196]; ABeh2=x[197]; ABel2=x[198] #AB2
  IBch2=x[199]; IBcl2=x[200]; IBah2=x[201]; IBal2=x[202]; IBeh2=x[203]; IBel2=x[204] #IB2
  RBch2=x[205]; RBcl2=x[206]; RBah2=x[207]; RBal2=x[208]; RBeh2=x[209]; RBel2=x[210] #RB2
  DBch2=x[211]; DBcl2=x[212]; DBah2=x[213]; DBal2=x[214]; DBeh2=x[215]; DBel2=x[216] #DB2
  
  
  with(as.list(params),{
    
    e_ch=betah
    e_cl=betal
    e_ah=betah 
    e_al=betal
    e_eh=betah
    e_el=betal
    e_vch=veic*betah  
    e_vcl=veic*betal 
    e_vah=vei*betah 
    e_val=vei*betal
    e_veh=vei*betah 
    e_vel=vei*betal 
    e_bch=vei2*betah  
    e_bcl=vei2*betal 
    e_bah=vei2*betah 
    e_bal=vei2*betal
    e_beh=vei2*betah 
    e_bel=vei2*betal
    
    #Update with CM by country
    #Malaysia
    #CM<-matrix(data=c(11.98, 3.77, 1.30, 6.14, 11.36, 2.03, 0.146, 0.141, 0.396), nrow=3)
    
    #India
    #CM<-matrix(data=c(15.78, 7.54, 1.61, 10.17, 24.73, 3.87, 0.191, 0.339, 0.254), nrow=3)
    
    #Ecuador
    CM<-matrix(data=c(11.98, 3.77, 1.3, 6.14, 11.36, 2.03, 0.146, 0.141, 0.396), nrow=3)
    #Sum of compartments, leaving D out
    Nch=sum(Sch, Ech, Ach, Ich, Rch)
    Nch2=sum(Sch2, Ech2, Ach2, Ich2, Rch2)
    NVch=sum(SVch, EVch, AVch, IVch, RVch)
    NVch2=sum(SVch2, EVch2, AVch2, IVch2, RVch2)
    NBch=sum(SBch, EBch, ABch, IBch, RBch)
    NBch2=sum(SBch2, EBch2, ABch2, IBch2, RBch2)
    
    Ncl=sum(Scl, Ecl, Acl, Icl, Rcl)
    Ncl2=sum(Scl2, Ecl2, Acl2, Icl2, Rcl2)
    NVcl=sum(SVcl, EVcl, AVcl, IVcl, RVcl)
    NVcl2=sum(SVcl2, EVcl2, AVcl2, IVcl2, RVcl2)
    NBcl=sum(SBcl, EBcl, ABcl, IBcl, RBcl)
    NBcl2=sum(SBcl2, EBcl2, ABcl2, IBcl2, RBcl2)
    
    Nah=sum(Sah, Eah, Aah, Iah, Rah)
    Nah2=sum(Sah2, Eah2, Aah2, Iah2, Rah2)
    NVah=sum(SVah, EVah, AVah, IVah, RVah)
    NVah2=sum(SVah2, EVah2, AVah2, IVah2, RVah2)
    NBah=sum(SBah, EBah, ABah, IBah, RBah)
    NBah2=sum(SBah2, EBah2, ABah2, IBah2, RBah2)
    
    Nal=sum(Sal, Eal, Aal, Ial, Ral)
    Nal2=sum(Sal2, Eal2, Aal2, Ial2, Ral2)
    NVal=sum(SVal, EVal, AVal, IVal, RVal)
    NVal2=sum(SVal2, EVal2, AVal2, IVal2, RVal2)
    NBal=sum(SBal, EBal, ABal, IBal, RBal)
    NBal2=sum(SBal2, EBal2, ABal2, IBal2, RBal2)
    
    Neh=sum(Seh, Eeh, Aeh, Ieh, Reh)
    Neh2=sum(Seh2, Eeh2, Aeh2, Ieh2, Reh2)
    NVeh=sum(SVeh, EVeh, AVeh, IVeh, RVeh)
    NVeh2=sum(SVeh2, EVeh2, AVeh2, IVeh2, RVeh2)
    NBeh=sum(SBeh, EBeh, ABeh, IBeh, RBeh)
    NBeh2=sum(SBeh2, EBeh2, ABeh2, IBeh2, RBeh2)
    
    Nel=sum(Sel, Eel, Ael, Iel, Rel)
    Nel2=sum(Sel2, Eel2, Ael2, Iel2, Rel2)
    NVel=sum(SVel, EVel, AVel, IVel, RVel)
    NVel2=sum(SVel2, EVel2, AVel2, IVel2, RVel2)
    NBel=sum(SBel, EBel, ABel, IBel, RBel)
    NBel2=sum(SBel2, EBel2, ABel2, IBel2, RBel2)
    
    Nchild=sum(Nch, Nch2, NVch, NVch2, NBch, NBch2, Ncl, Ncl2, NVcl, NVcl2, NBcl, NBcl2)
    Nadult=sum(Nah, Nah2, NVah, NVah2, NBah, NBah2, Nal, Nal2, NVal, NVal2, NBal, NBal2)
    Nold=sum(Neh, Neh2, NVeh, NVeh2, NBeh, NBeh2, Nel, Nel2, NVel, NVel2, NBel, NBel2)
    
    #Avoid division by zero, must do after summation for age stratified Ns
    if(Nch==0){Nch<-0.01}
    if(Ncl==0){Ncl<-0.01}
    if(Nah==0){Nah<-0.01}
    if(Nal==0){Nal<-0.01}
    if(Neh==0){Neh<-0.01}
    if(Nel==0){Nel<-0.01}
    if(Nch2==0){Nch2<-0.01}
    if(Ncl2==0){Ncl2<-0.01}
    if(Nah2==0){Nah2<-0.01}
    if(Nal2==0){Nal2<-0.01}
    if(Neh2==0){Neh2<-0.01}
    if(Nel2==0){Nel2<-0.01}
    if(NVch==0){NVch<-0.01}
    if(NVcl==0){NVcl<-0.01}
    if(NVah==0){NVah<-0.01}
    if(NVal==0){NVal<-0.01}
    if(NVeh==0){NVeh<-0.01}
    if(NVel==0){NVel<-0.01}
    if(NVch2==0){NVch2<-0.01}
    if(NVcl2==0){NVcl2<-0.01}
    if(NVah2==0){NVah2<-0.01}
    if(NVal2==0){NVal2<-0.01}
    if(NVeh2==0){NVeh2<-0.01}
    if(NVel2==0){NVel2<-0.01}
    if(NBch==0){NBch<-0.01}
    if(NBcl==0){NBcl<-0.01}
    if(NBah==0){NBah<-0.01}
    if(NBal==0){NBal<-0.01}
    if(NBeh==0){NBeh<-0.01}
    if(NBel==0){NBel<-0.01}
    if(NBch2==0){NBch2<-0.01}
    if(NBcl2==0){NBcl2<-0.01}
    if(NBah2==0){NBah2<-0.01}
    if(NBal2==0){NBal2<-0.01}
    if(NBeh2==0){NBeh2<-0.01}
    if(NBel2==0){NBel2<-0.01}
    
    #infect.child <- Ich+Icl+Ach+Acl+Ich2+Icl2+Ach2+Acl2+IVch+IVcl+AVch+AVcl
    #infect.adult <- Iah+Ial+Aah+Aal+Iah2+Ial2+Aah2+Aal2+IVah+IVal+AVah+AVal
    #infect.old   <- Ieh+Iel+Aeh+Ael+Ieh2+Iel2+Aeh2+Ael2+IVeh+IVel+AVeh+AVel
    
    scale_immune <- 0.27
    
    
    infect.child <- Ich+Icl+Ach+Acl+ ## not infected or vaccinated
      scale_immune*(Ich2+Icl2+Ach2+Acl2+ ## infected once, not yet vaccinated
      IVch+IVcl+AVch+AVcl + ## vaccinated, not infected
      IVch2+IVcl2+AVch2+AVcl2 + ## vaccinated, infected
      IBch + IBcl + ABch + ABcl + ## boosted, not infected
      IBch2 + IBcl2 + ABch2 + ABcl2) ## boosted, infected
    
    infect.adult <- Iah+Ial+Aah+Aal+ ## not infected or vaccinated
      scale_immune*(Iah2+Ial2+Aah2+Aal2+ ## infected once, not yet vaccinated
      IVah+IVal+AVah+AVal + ## vaccinated, not infected
      IVah2+IVal2+AVah2+AVal2 + ## vaccinated, infected
      IBah + IBal + ABah + ABal + ## boosted, not infected
      IBah2 + IBal2 + ABah2 + ABal2) ## boosted, infected
    
    infect.old   <- Ieh+Iel+Aeh+Ael+ ## not infected or vaccinated
      scale_immune*(Ieh2+Iel2+Aeh2+Ael2+ ## infected once, not yet vaccinated
      IVeh+IVel+AVeh+AVel + ## vaccinated, not infected
      IVeh2+IVel2+AVeh2+AVel2 + ## vaccinated, infected
      IBeh + IBel + ABeh + ABel + ## boosted, not infected
      IBeh2 + IBel2 + ABeh2 + ABel2) ## boosted, infected
    
    
    #Vaccination is determined by rate
    #Primary vaccination, assumed to be zero during this wave (see code below to change if desired)
    lambda_Sch<-0
    lambda_Scl<-0
    lambda_Reh<-0; lambda_Seh<-0; lambda_Rel<-0; lambda_Sel<-0
    lambda_Rah<-0; lambda_Sah<-0; lambda_Ral<-0; lambda_Sal<-0
    lambda_Rch<-lambda_Sch; lambda_Rcl<-lambda_Scl
    
    #Booster vaccination.  For these runs, we have no child boosting (boost.child=0 across all parameter sets). 
    #Keeping code available so we can change later if needed.
    
    ## notes: lambda_true * (Sv + Ev + Rv) = lambda_calcd * (Sv + Ev + Iv + Av + Rv + Sb + Eb + Ib + Ab + Rb)
    ## lambda_true = lambda_calcd * (Sv + Ev + Iv + Av + Rv + Sb + Eb + Ib + Ab + Rb)/(Sv + Ev + Rv)
    ## Sv + Ev + Iv + Av + Rv + Sb + Eb + Ib + Ab + Rb is the total number boosted
    ## Sv + Ev + Rv is the number eligible to be vaccinated
    ## will need to set a cieling: min(1, lambda_true)
    
    if(boost.child==0){
      lambda_Sch2<-0
      lambda_Scl2<-0
    }
    if(boost.child>0){
      lambda_Sch2<-min(c(vacc_H(t+boost.start+1)*(NVch + NVch2 + NBch + NBch2)/(NVch+NVch2-(IVch+IVch2+AVch+AVch2)),1))
      lambda_Scl2<-min(c(vacc_L(t+boost.start+1)*(NVcl + NVcl2 + NBcl + NBcl2)/(NVcl+NVcl2-(IVcl+IVcl2+AVcl+AVcl2)),1))
    }
    if(boost.adult>0){
      lambda_Sah2<-min(c(vacc_H(t+boost.start+1) * (NVah + NVah2 + NBah + NBah2)/(NVah+NVah2-(IVah+IVah2+AVah+AVah2)),1))
      lambda_Sal2<-min(c(vacc_L(t+boost.start+1) * (NVal + NVal2 + NBal + NBal2)/(NVal+NVal2-(IVal+IVal2+AVal+AVal2)),1))
      lambda_Seh2<-min(c(vacc_H(t+boost.start+1) * (NVeh + NVeh2 + NBeh + NBeh2)/(NVeh+NVeh2-(IVeh+IVeh2+AVeh+AVeh2)),1))
      lambda_Sel2<-min(c(vacc_L(t+boost.start+1) * (NVel + NVel2 + NBel + NBel2)/(NVel+NVel2-(IVel+IVel2+AVel+AVel2)),1))
    }
    if(boost.adult==0){
      lambda_Sah2<-0
      lambda_Sal2<-0
      lambda_Seh2<-0
      lambda_Sel2<-0
    }
    
    #Old code for child primary vaccination, would need to account for rollout rates using type III functions if want to add back
    #This code may be helpful in that case so keeping it but commenting out
    #if(cvax>0){
    #  if(((NVch+NVch2+NBch+NBch2)/(Nch+Nch2+NVch+NVch2+NBch+NBch2))>child.max){ #If coverage has peaked in a group, stop vaccinating
    #    lambda_Sch<-0}
    #  if(((NVch+NVch2+NBch+NBch2)/(Nch+Nch2+NVch+NVch2+NBch+NBch2))<=child.max){ #If coverage has peaked in a group, allocate according to weights
    #    lambda_Sch<-childrate.h}
      #Same thing for cl
    #  if(((NVcl+NVcl2+NBcl+NBcl2)/(Ncl+Ncl2+NVcl+NVcl2+NBcl+NBcl2))>child.max){
    #    lambda_Scl<-0}
    #  if(((NVcl+NVcl2+NBcl+NBcl2)/(Ncl+Ncl2+NVcl+NVcl2+NBcl+NBcl2))<=child.max){
    #    lambda_Scl<-childrate.l}      
    #}
    
    #if(cvax==0){
    #  lambda_Sch<-0
    #  lambda_Scl<-0}
    
    #if(boost.child>0){
    #  if(((NBch+NBcl2)/(NVch+NVch2+NBch+NBch2))>max.boost.child){ 
    #    lambda_Sch2<-0} 
    #  if(((NBch+NBch2)/(NVch+NVch2+NBch+NBch2))<=max.boost.child){ #If coverage is less than 80%, allocate according to weights
    #    lambda_Sch2<-vocrate.h}
    #  if(((NBcl+NBcl2)/(NVcl+NVcl2+NBcl+NBcl2))>max.boost.child){
    #    lambda_Scl2<-0}
    #  if(((NBcl+NBcl2)/(NVcl+NVcl2+NBcl+NBcl2))<=max.boost.child){
    #    lambda_Scl2<-vocrate.l}}
    
    
    #if(boost.child==0){
    #  lambda_Sch2<-0
    #  lambda_Scl2<-0
    #}
    
    lambda_Rch2<-lambda_Sch2; lambda_Rcl2<-lambda_Scl2
    lambda_Rah2<-lambda_Sah2; lambda_Ral2<-lambda_Sal2
    lambda_Reh2<-lambda_Seh2; lambda_Rel2<-lambda_Sel2
    
    contactC <- (CM[1,1]/Nchild*(infect.child)) + 
      (CM[1,2]/Nadult*(infect.adult)) + (CM[1,3]/Nold*(infect.old))
    
    contactA <- (CM[2,1]/Nchild*(infect.child)) + 
      (CM[2,2]/Nadult*(infect.adult))+ (CM[2,3]/Nold*(infect.old))
    
    contactO <- (CM[3,1]/Nchild*(infect.child)) + 
      (CM[3,2]/Nadult*(infect.adult)) + (CM[3,3]/Nold*(infect.old))
    
    ################################## 
    ################################## 
    
    ################################## 
    #Never infected or vaccinated
    #Susceptible
    dSch = -Sch*e_ch*sd*contactC - lambda_Sch*Sch
    dScl = -Scl*e_cl*sd*contactC - lambda_Scl*Scl
    dSah = -Sah*e_ah*sd*contactA - lambda_Sah*Sah
    dSal = -Sal*e_al*sd*contactA - lambda_Sal*Sal
    dSeh = -Seh*e_eh*sd*contactO - lambda_Seh*Seh
    dSel = -Sel*e_el*sd*contactO - lambda_Sel*Sel
    
    #Exposed 
    dEch = Sch*e_ch*sd*contactC - sigma_Ech*Ech - lambda_Sch*Ech
    dEcl = Scl*e_cl*sd*contactC - sigma_Ecl*Ecl - lambda_Scl*Ecl
    dEah = Sah*e_ah*sd*contactA - sigma_Eah*Eah - lambda_Sah*Eah
    dEal = Sal*e_al*sd*contactA - sigma_Eal*Eal - lambda_Sal*Eal
    dEeh = Seh*e_eh*sd*contactO - sigma_Eeh*Eeh - lambda_Seh*Eeh
    dEel = Sel*e_el*sd*contactO - sigma_Eel*Eel - lambda_Sel*Eel
    
    #Asymptomatic                                     
    dAch = (1-nu_Ech)*sigma_Ech*Ech -gA_Ach*Ach
    dAcl = (1-nu_Ecl)*sigma_Ecl*Ecl -gA_Acl*Acl
    dAah = (1-nu_Eah)*sigma_Eah*Eah -gA_Aah*Aah
    dAal = (1-nu_Eal)*sigma_Eal*Eal -gA_Aal*Aal
    dAeh = (1-nu_Eeh)*sigma_Eeh*Eeh -gA_Aeh*Aeh
    dAel = (1-nu_Eel)*sigma_Eel*Eel -gA_Ael*Ael
    
    #Symptomatic
    dIch = nu_Ech*sigma_Ech*Ech -gI_Ich*Ich
    dIcl = nu_Ecl*sigma_Ecl*Ecl -gI_Icl*Icl
    dIah = nu_Eah*sigma_Eah*Eah -gI_Iah*Iah
    dIal = nu_Eal*sigma_Eal*Eal -gI_Ial*Ial
    dIeh = nu_Eeh*sigma_Eeh*Eeh -gI_Ieh*Ieh
    dIel = nu_Eel*sigma_Eel*Eel -gI_Iel*Iel
    
    #Recovered 
    dRch = gA_Ach*Ach + (1-rho_Ich)*gI_Ich*Ich - lambda_Rch*Rch - omega*Rch
    dRcl = gA_Acl*Acl + (1-rho_Icl)*gI_Icl*Icl - lambda_Rcl*Rcl - omega*Rcl 
    dRah = gA_Aah*Aah + (1-rho_Iah)*gI_Iah*Iah - lambda_Rah*Rah - omega*Rah
    dRal = gA_Aal*Aal + (1-rho_Ial)*gI_Ial*Ial - lambda_Ral*Ral - omega*Ral
    dReh = gA_Aeh*Aeh + (1-rho_Ieh)*gI_Ieh*Ieh - lambda_Reh*Reh - omega*Reh
    dRel = gA_Ael*Ael + (1-rho_Iel)*gI_Iel*Iel - lambda_Rel*Rel - omega*Rel
    
    dDch = rho_Ich*gI_Ich*Ich
    dDcl = rho_Icl*gI_Icl*Icl
    dDah = rho_Iah*gI_Iah*Iah
    dDal = rho_Ial*gI_Ial*Ial
    dDeh = rho_Ieh*gI_Ieh*Ieh
    dDel = rho_Iel*gI_Iel*Iel
    
    ################################## 
    #Infected once, not yet vaccinated
    #Susceptible
    dSch2 = -Sch2*e_ch*sd*(1-(0.8*cp))*contactC - lambda_Sch*Sch2 + omega*(Rch+Rch2)
    dScl2 = -Scl2*e_cl*sd*(1-(0.8*cp))*contactC - lambda_Scl*Scl2 + omega*(Rcl+Rcl2)
    dSah2 = -Sah2*e_ah*sd*(1-(0.8*cp))*contactA - lambda_Sah*Sah2 + omega*(Rah+Rah2)
    dSal2 = -Sal2*e_al*sd*(1-(0.8*cp))*contactA - lambda_Sal*Sal2 + omega*(Ral+Ral2)
    dSeh2 = -Seh2*e_eh*sd*(1-(0.8*cp))*contactO - lambda_Seh*Seh2 + omega*(Reh+Reh2)
    dSel2 = -Sel2*e_el*sd*(1-(0.8*cp))*contactO - lambda_Sel*Sel2 + omega*(Rel+Rel2)
    
    #Exposed 
    dEch2 = Sch2*e_ch*sd*(1-(0.8*cp))*contactC - sigma_Ech*Ech2 - lambda_Sch*Ech2
    dEcl2 = Scl2*e_cl*sd*(1-(0.8*cp))*contactC - sigma_Ecl*Ecl2 - lambda_Scl*Ecl2
    dEah2 = Sah2*e_ah*sd*(1-(0.8*cp))*contactA - sigma_Eah*Eah2 - lambda_Sah*Eah2
    dEal2 = Sal2*e_al*sd*(1-(0.8*cp))*contactA - sigma_Eal*Eal2 - lambda_Sal*Eal2
    dEeh2 = Seh2*e_eh*sd*(1-(0.8*cp))*contactO - sigma_Eeh*Eeh2 - lambda_Seh*Eeh2
    dEel2 = Sel2*e_el*sd*(1-(0.8*cp))*contactO - sigma_Eel*Eel2 - lambda_Sel*Eel2
    
    #Asymptomatic                                     
    dAch2 = (1-nu_Ech)*sigma_Ech*Ech2 - gA_Ach*Ach2
    dAcl2 = (1-nu_Ecl)*sigma_Ecl*Ecl2 - gA_Acl*Acl2
    dAah2 = (1-nu_Eah)*sigma_Eah*Eah2 - gA_Aah*Aah2
    dAal2 = (1-nu_Eal)*sigma_Eal*Eal2 - gA_Aal*Aal2
    dAeh2 = (1-nu_Eeh)*sigma_Eeh*Eeh2 - gA_Aeh*Aeh2
    dAel2 = (1-nu_Eel)*sigma_Eel*Eel2 - gA_Ael*Ael2
    
    #Symptomatic
    dIch2 = nu_Ech*sigma_Ech*Ech2 - gI_Ich*Ich2
    dIcl2 = nu_Ecl*sigma_Ecl*Ecl2 - gI_Icl*Icl2
    dIah2 = nu_Eah*sigma_Eah*Eah2 - gI_Iah*Iah2
    dIal2 = nu_Eal*sigma_Eal*Eal2 - gI_Ial*Ial2
    dIeh2 = nu_Eeh*sigma_Eeh*Eeh2 - gI_Ieh*Ieh2
    dIel2 = nu_Eel*sigma_Eel*Eel2 - gI_Iel*Iel2
    
    #Recovered 
    dRch2 = gA_Ach*Ach2 + (1-((1-cp)*rho_Ich))*gI_Ich*Ich2 - lambda_Rch*Rch2 - omega*Rch2
    dRcl2 = gA_Acl*Acl2 + (1-((1-cp)*rho_Icl))*gI_Icl*Icl2 - lambda_Rcl*Rcl2 - omega*Rcl2 
    dRah2 = gA_Aah*Aah2 + (1-((1-cp)*rho_Iah))*gI_Iah*Iah2 - lambda_Rah*Rah2 - omega*Rah2
    dRal2 = gA_Aal*Aal2 + (1-((1-cp)*rho_Ial))*gI_Ial*Ial2 - lambda_Ral*Ral2 - omega*Ral2
    dReh2 = gA_Aeh*Aeh2 + (1-((1-cp)*rho_Ieh))*gI_Ieh*Ieh2 - lambda_Reh*Reh2 - omega*Reh2
    dRel2 = gA_Ael*Ael2 + (1-((1-cp)*rho_Iel))*gI_Iel*Iel2 - lambda_Rel*Rel2 - omega*Rel2
    
    dDch2 = rho_Ich*(1-cp)*gI_Ich*Ich2
    dDcl2 = rho_Icl*(1-cp)*gI_Icl*Icl2
    dDah2 = rho_Iah*(1-cp)*gI_Iah*Iah2
    dDal2 = rho_Ial*(1-cp)*gI_Ial*Ial2
    dDeh2 = rho_Ieh*(1-cp)*gI_Ieh*Ieh2
    dDel2 = rho_Iel*(1-cp)*gI_Iel*Iel2
    
    ##################################
    #Uninfected, vaccinated   
    #Susceptible
    dSVch = -SVch*e_vch*sd*contactC + lambda_Sch*Sch - lambda_Sch2*SVch 
    dSVcl = -SVcl*e_vcl*sd*contactC + lambda_Scl*Scl - lambda_Scl2*SVcl
    dSVah = -SVah*e_vah*sd*contactA + lambda_Sah*Sah - lambda_Sah2*SVah
    dSVal = -SVal*e_val*sd*contactA + lambda_Sal*Sal - lambda_Sal2*SVal
    dSVeh = -SVeh*e_veh*sd*contactO + lambda_Seh*Seh - lambda_Seh2*SVeh
    dSVel = -SVel*e_vel*sd*contactO + lambda_Sel*Sel - lambda_Sel2*SVel
    
    #Exposed, vaccinated
    dEVch = SVch*e_vch*sd*contactC - sigmav_EVch*EVch + lambda_Sch*Ech - lambda_Sch2*EVch
    dEVcl = SVcl*e_vcl*sd*contactC - sigmav_EVcl*EVcl + lambda_Scl*Ecl - lambda_Scl2*EVcl
    dEVah = SVah*e_vah*sd*contactA - sigmav_EVah*EVah + lambda_Sah*Eah - lambda_Sah2*EVah
    dEVal = SVal*e_val*sd*contactA - sigmav_EVal*EVal + lambda_Sal*Eal - lambda_Sal2*EVal
    dEVeh = SVeh*e_veh*sd*contactO - sigmav_EVeh*EVeh + lambda_Seh*Eeh - lambda_Seh2*EVeh
    dEVel = SVel*e_vel*sd*contactO - sigmav_EVel*EVel + lambda_Sel*Eel - lambda_Sel2*EVel
    
    #Asymptomatic, vaccinated                                                      
    dAVch = (1-nuv_EVch)*sigmav_EVch*EVch - gAV_AVch*AVch
    dAVcl = (1-nuv_EVcl)*sigmav_EVcl*EVcl - gAV_AVcl*AVcl
    dAVah = (1-nuv_EVah)*sigmav_EVah*EVah - gAV_AVah*AVah
    dAVal = (1-nuv_EVal)*sigmav_EVal*EVal - gAV_AVal*AVal
    dAVeh = (1-nuv_EVeh)*sigmav_EVeh*EVeh - gAV_AVeh*AVeh
    dAVel = (1-nuv_EVel)*sigmav_EVel*EVel - gAV_AVel*AVel
    
    #Symptomatic, vaccinated
    dIVch = nuv_EVch*sigmav_EVch*EVch - gIV_IVch*IVch
    dIVcl = nuv_EVcl*sigmav_EVcl*EVcl - gIV_IVcl*IVcl
    dIVah = nuv_EVah*sigmav_EVah*EVah - gIV_IVah*IVah
    dIVal = nuv_EVal*sigmav_EVal*EVal - gIV_IVal*IVal
    dIVeh = nuv_EVeh*sigmav_EVeh*EVeh - gIV_IVeh*IVeh
    dIVel = nuv_EVel*sigmav_EVel*EVel - gIV_IVel*IVel
    
    #Recovered, vaccinated
    dRVch = gAV_AVch*AVch + (1-rhov_IVch)*gIV_IVch*IVch + lambda_Rch*Rch - 
      omega*RVch - lambda_Rch2*RVch 
    dRVcl = gAV_AVcl*AVcl + (1-rhov_IVcl)*gIV_IVcl*IVcl + lambda_Rcl*Rcl - 
      omega*RVcl - lambda_Rcl2*RVcl 
    dRVah = gAV_AVah*AVah + (1-rhov_IVah)*gIV_IVah*IVah + lambda_Rah*Rah - 
      omega*RVah - lambda_Rah2*RVah 
    dRVal = gAV_AVal*AVal + (1-rhov_IVal)*gIV_IVal*IVal + lambda_Ral*Ral - 
      omega*RVal - lambda_Ral2*RVal 
    dRVeh = gAV_AVeh*AVeh + (1-rhov_IVeh)*gIV_IVeh*IVeh + lambda_Reh*Reh - 
      omega*RVeh - lambda_Reh2*RVeh 
    dRVel = gAV_AVel*AVel + (1-rhov_IVel)*gIV_IVel*IVel + lambda_Rel*Rel - 
      omega*RVel - lambda_Rel2*RVel 
    
    dDVch = rhov_IVch*gIV_IVch*IVch
    dDVcl = rhov_IVcl*gIV_IVcl*IVcl
    dDVah = rhov_IVah*gIV_IVah*IVah
    dDVal = rhov_IVal*gIV_IVal*IVal
    dDVeh = rhov_IVeh*gIV_IVeh*IVeh
    dDVel = rhov_IVel*gIV_IVel*IVel
    
    ################################## 
    #Primary vaccination + 1 natural infection
    #Susceptible
    dSVch2 = -SVch2*e_vch*sd*(1-(0.8*cp))*contactC + lambda_Sch*Sch2 + omega*(RVch+RVch2) - lambda_Sch2*SVch2 
    dSVcl2 = -SVcl2*e_vcl*sd*(1-(0.8*cp))*contactC + lambda_Scl*Scl2 + omega*(RVcl+RVcl2) - lambda_Scl2*SVcl2 
    dSVah2 = -SVah2*e_vah*sd*(1-(0.8*cp))*contactA + lambda_Sah*Sah2 + omega*(RVah+RVah2) - lambda_Sah2*SVah2 
    dSVal2 = -SVal2*e_val*sd*(1-(0.8*cp))*contactA + lambda_Sal*Sal2 + omega*(RVal+RVal2) - lambda_Sal2*SVal2 
    dSVeh2 = -SVeh2*e_veh*sd*(1-(0.8*cp))*contactO + lambda_Seh*Seh2 + omega*(RVeh+RVeh2) - lambda_Seh2*SVeh2 
    dSVel2 = -SVel2*e_vel*sd*(1-(0.8*cp))*contactO + lambda_Sel*Sel2 + omega*(RVel+RVel2) - lambda_Sel2*SVel2 
    
    #Exposed
    dEVch2 = SVch2*e_vch*sd*(1-(0.8*cp))*contactC - sigmav_EVch*EVch2 + lambda_Sch*Ech2 - lambda_Sch2*EVch2
    dEVcl2 = SVcl2*e_vcl*sd*(1-(0.8*cp))*contactC - sigmav_EVcl*EVcl2 + lambda_Scl*Ecl2 - lambda_Scl2*EVcl2
    dEVah2 = SVah2*e_vah*sd*(1-(0.8*cp))*contactA - sigmav_EVah*EVah2 + lambda_Sah*Eah2 - lambda_Sah2*EVah2
    dEVal2 = SVal2*e_val*sd*(1-(0.8*cp))*contactA - sigmav_EVal*EVal2 + lambda_Sal*Eal2 - lambda_Sal2*EVal2
    dEVeh2 = SVeh2*e_veh*sd*(1-(0.8*cp))*contactO - sigmav_EVeh*EVeh2 + lambda_Seh*Eeh2 - lambda_Seh2*EVeh2
    dEVel2 = SVel2*e_vel*sd*(1-(0.8*cp))*contactO - sigmav_EVel*EVel2 + lambda_Sel*Eel2 - lambda_Sel2*EVel2
    
    #Asymptomatic                                                      
    dAVch2 = (1-nuv_EVch)*sigmav_EVch*EVch2 - gAV_AVch*AVch2
    dAVcl2 = (1-nuv_EVcl)*sigmav_EVcl*EVcl2 - gAV_AVcl*AVcl2
    dAVah2 = (1-nuv_EVah)*sigmav_EVah*EVah2 - gAV_AVah*AVah2
    dAVal2 = (1-nuv_EVal)*sigmav_EVal*EVal2 - gAV_AVal*AVal2
    dAVeh2 = (1-nuv_EVeh)*sigmav_EVeh*EVeh2 - gAV_AVeh*AVeh2
    dAVel2 = (1-nuv_EVel)*sigmav_EVel*EVel2 - gAV_AVel*AVel2
    
    #Symptomatic, vaccinated
    dIVch2 = nuv_EVch*sigmav_EVch*EVch2 - gIV_IVch*IVch2
    dIVcl2 = nuv_EVcl*sigmav_EVcl*EVcl2 - gIV_IVcl*IVcl2
    dIVah2 = nuv_EVah*sigmav_EVah*EVah2 - gIV_IVah*IVah2
    dIVal2 = nuv_EVal*sigmav_EVal*EVal2 - gIV_IVal*IVal2
    dIVeh2 = nuv_EVeh*sigmav_EVeh*EVeh2 - gIV_IVeh*IVeh2
    dIVel2 = nuv_EVel*sigmav_EVel*EVel2 - gIV_IVel*IVel2
    
    #Recovered, vaccinated
    dRVch2 = gAV_AVch*AVch2 + (1-((1-cp)*rhov_IVch))*gIV_IVch*IVch2 - 
      omega*RVch2 + lambda_Rch*Rch2 - lambda_Rch2*RVch2 
    dRVcl2 = gAV_AVcl*AVcl2 + (1-((1-cp)*rhov_IVcl))*gIV_IVcl*IVcl2 - 
      omega*RVcl2 + lambda_Rcl*Rcl2 - lambda_Rcl2*RVcl2
    dRVah2 = gAV_AVah*AVah2 + (1-((1-cp)*rhov_IVah))*gIV_IVah*IVah2 - 
      omega*RVah2 + lambda_Rah*Rah2 - lambda_Rah2*RVah2
    dRVal2 = gAV_AVal*AVal2 + (1-((1-cp)*rhov_IVal))*gIV_IVal*IVal2 - 
      omega*RVal2 + lambda_Ral*Ral2 - lambda_Ral2*RVal2
    dRVeh2=  gAV_AVeh*AVeh2 + (1-((1-cp)*rhov_IVeh))*gIV_IVeh*IVeh2 - 
      omega*RVeh2 + lambda_Reh*Reh2 - lambda_Reh2*RVeh2
    dRVel2 = gAV_AVel*AVel2 + (1-((1-cp)*rhov_IVel))*gIV_IVel*IVel2 - 
      omega*RVel2 + lambda_Rel*Rel2 - lambda_Rel2*RVel2
    
    dDVch2 = (1-cp)*rhov_IVch*gIV_IVch*IVch2
    dDVcl2 = (1-cp)*rhov_IVcl*gIV_IVcl*IVcl2
    dDVah2 = (1-cp)*rhov_IVah*gIV_IVah*IVah2
    dDVal2 = (1-cp)*rhov_IVal*gIV_IVal*IVal2
    dDVeh2 = (1-cp)*rhov_IVeh*gIV_IVeh*IVeh2
    dDVel2 = (1-cp)*rhov_IVel*gIV_IVel*IVel2
    
    ################################## 
    #Boosted and uninfected
    #Susceptible
    dSBch = -SBch*e_bch*sd*contactC + lambda_Sch2*SVch
    dSBcl = -SBcl*e_bcl*sd*contactC + lambda_Scl2*SVcl
    dSBah = -SBah*e_bah*sd*contactA + lambda_Sah2*SVah
    dSBal = -SBal*e_bal*sd*contactA + lambda_Sal2*SVal
    dSBeh = -SBeh*e_beh*sd*contactO + lambda_Seh2*SVeh
    dSBel = -SBel*e_bel*sd*contactO + lambda_Sel2*SVel
    
    #Exposed, boosted
    dEBch = SBch*e_bch*sd*contactC - sigmav_EVch*EBch + lambda_Sch2*EVch
    dEBcl = SBcl*e_bcl*sd*contactC - sigmav_EVcl*EBcl + lambda_Scl2*EVcl
    dEBah = SBah*e_bah*sd*contactA - sigmav_EVah*EBah + lambda_Sah2*EVah
    dEBal = SBal*e_bal*sd*contactA - sigmav_EVal*EBal + lambda_Sal2*EVal
    dEBeh = SBeh*e_beh*sd*contactO - sigmav_EVeh*EBeh + lambda_Seh2*EVeh
    dEBel = SBel*e_bel*sd*contactO - sigmav_EVel*EBel + lambda_Sel2*EVel
    
    #Asymptomatic, boosted
    dABch = (1-nuv_EVch)*sigmav_EVch*EBch - gAV_AVch*ABch
    dABcl = (1-nuv_EVcl)*sigmav_EVcl*EBcl - gAV_AVcl*ABcl
    dABah = (1-nuv_EVah)*sigmav_EVah*EBah - gAV_AVah*ABah
    dABal = (1-nuv_EVal)*sigmav_EVal*EBal - gAV_AVal*ABal
    dABeh = (1-nuv_EVeh)*sigmav_EVeh*EBeh - gAV_AVeh*ABeh
    dABel = (1-nuv_EVel)*sigmav_EVel*EBel - gAV_AVel*ABel
    
    #Symptomatic, boosted
    dIBch = nuv_EVch*sigmav_EVch*EBch - gIV_IVch*IBch 
    dIBcl = nuv_EVcl*sigmav_EVcl*EBcl - gIV_IVcl*IBcl 
    dIBah = nuv_EVah*sigmav_EVah*EBah - gIV_IVah*IBah 
    dIBal = nuv_EVal*sigmav_EVal*EBal - gIV_IVal*IBal 
    dIBeh = nuv_EVeh*sigmav_EVeh*EBeh - gIV_IVeh*IBeh 
    dIBel = nuv_EVel*sigmav_EVel*EBel - gIV_IVel*IBel 
    
    #Recovered, boosted
    dRBch = gAV_AVch*ABch + (1-rhov_IVch2)*gIV_IVch*IBch + 
      lambda_Rch2*RVch - omega*RBch 
    dRBcl = gAV_AVcl*ABcl + (1-rhov_IVcl2)*gIV_IVcl*IBcl + 
      lambda_Rcl2*RVcl - omega*RBcl 
    dRBah = gAV_AVah*ABah + (1-rhov_IVah2)*gIV_IVah*IBah + 
      lambda_Rah2*RVah - omega*RBah 
    dRBal = gAV_AVal*ABal + (1-rhov_IVal2)*gIV_IVal*IBal + 
      lambda_Ral2*RVal - omega*RBal 
    dRBeh = gAV_AVeh*ABeh + (1-rhov_IVeh2)*gIV_IVeh*IBeh + 
      lambda_Reh2*RVeh - omega*RBeh 
    dRBel = gAV_AVel*ABel + (1-rhov_IVel2)*gIV_IVel*IBel + 
      lambda_Rel2*RVel - omega*RBel 
    
    #Deaths, boosted
    dDBch = rhov_IVch2*gIV_IVch*IBch
    dDBcl = rhov_IVcl2*gIV_IVcl*IBcl
    dDBah = rhov_IVah2*gIV_IVah*IBah
    dDBal = rhov_IVal2*gIV_IVal*IBal
    dDBeh = rhov_IVeh2*gIV_IVeh*IBeh
    dDBel = rhov_IVel2*gIV_IVel*IBel
    
    ################################## 
    #Susceptible, boosted, 1 natural infection
    #Susceptible
    dSBch2 = -SBch2*e_bch*sd*(1-(0.8*cp))*contactC + lambda_Sch2*SVch2 + omega*(RBch+RBch2) 
    dSBcl2 = -SBcl2*e_bcl*sd*(1-(0.8*cp))*contactC + lambda_Scl2*SVcl2 + omega*(RBcl+RBcl2) 
    dSBah2 = -SBah2*e_bah*sd*(1-(0.8*cp))*contactA + lambda_Sah2*SVah2 + omega*(RBah+RBah2)
    dSBal2 = -SBal2*e_bal*sd*(1-(0.8*cp))*contactA + lambda_Sal2*SVal2 + omega*(RBal+RBal2)
    dSBeh2 = -SBeh2*e_beh*sd*(1-(0.8*cp))*contactO + lambda_Seh2*SVeh2 + omega*(RBeh+RBeh2) 
    dSBel2 = -SBel2*e_bel*sd*(1-(0.8*cp))*contactO + lambda_Sel2*SVel2 + omega*(RBel+RBel2)
    
    #Exposed, boosted
    dEBch2 = SBch2*e_bch*sd*(1-(0.8*cp))*contactC - sigmav_EVch*EBch2 + lambda_Sch2*EVch2
    dEBcl2 = SBcl2*e_bcl*sd*(1-(0.8*cp))*contactC - sigmav_EVcl*EBcl2 + lambda_Scl2*EVcl2
    dEBah2 = SBah2*e_bah*sd*(1-(0.8*cp))*contactA - sigmav_EVah*EBah2 + lambda_Sah2*EVah2 
    dEBal2 = SBal2*e_bal*sd*(1-(0.8*cp))*contactA - sigmav_EVal*EBal2 + lambda_Sal2*EVal2
    dEBeh2 = SBeh2*e_beh*sd*(1-(0.8*cp))*contactO - sigmav_EVeh*EBeh2 + lambda_Seh2*EVeh2
    dEBel2 = SBel2*e_bel*sd*(1-(0.8*cp))*contactO - sigmav_EVel*EBel2 + lambda_Sel2*EVel2
    
    #Asymptomatic, boosted
    dABch2 = (1-nuv_EVch)*sigmav_EVch*EBch2 - gAV_AVch*ABch2
    dABcl2 = (1-nuv_EVcl)*sigmav_EVcl*EBcl2 - gAV_AVcl*ABcl2
    dABah2 = (1-nuv_EVah)*sigmav_EVah*EBah2 - gAV_AVah*ABah2
    dABal2 = (1-nuv_EVal)*sigmav_EVal*EBal2 - gAV_AVal*ABal2
    dABeh2 = (1-nuv_EVeh)*sigmav_EVeh*EBeh2 - gAV_AVeh*ABeh2
    dABel2 = (1-nuv_EVel)*sigmav_EVel*EBel2 - gAV_AVel*ABel2
    
    #Symptomatic, boosted
    dIBch2 = nuv_EVch*sigmav_EVch*EBch2 - gIV_IVch*IBch2 
    dIBcl2 = nuv_EVcl*sigmav_EVcl*EBcl2 - gIV_IVcl*IBcl2 
    dIBah2 = nuv_EVah*sigmav_EVah*EBah2 - gIV_IVah*IBah2 
    dIBal2 = nuv_EVal*sigmav_EVal*EBal2 - gIV_IVal*IBal2 
    dIBeh2 = nuv_EVeh*sigmav_EVeh*EBeh2 - gIV_IVeh*IBeh2 
    dIBel2 = nuv_EVel*sigmav_EVel*EBel2 - gIV_IVel*IBel2 
    
    #Recovered, boosted
    dRBch2 = gAV_AVch*ABch2 +(1-(1-cp)*rhov_IVch2)*gIV_IVch*IBch2 + 
      lambda_Rch2*RVch2 - omega*RBch2 
    dRBcl2 = gAV_AVcl*ABcl2 +(1-(1-cp)*rhov_IVcl2)*gIV_IVcl*IBcl2 + 
      lambda_Rcl2*RVcl2 - omega*RBcl2 
    dRBah2 = gAV_AVah*ABah2 +(1-(1-cp)*rhov_IVah2)*gIV_IVah*IBah2 + 
      lambda_Rah2*RVah2 - omega*RBah2  
    dRBal2 = gAV_AVal*ABal2 +(1-(1-cp)*rhov_IVal2)*gIV_IVal*IBal2 + 
      lambda_Ral2*RVal2 - omega*RBal2  
    dRBeh2 = gAV_AVeh*ABeh2 +(1-(1-cp)*rhov_IVeh2)*gIV_IVeh*IBeh2 + 
      lambda_Reh2*RVeh2 - omega*RBeh2  
    dRBel2 = gAV_AVel*ABel2 +(1-(1-cp)*rhov_IVel2)*gIV_IVel*IBel2 + 
      lambda_Rel2*RVel2 - omega*RBel2  
    
    #Deaths, boosted
    dDBch2 = (1-cp)*rhov_IVch2*gIV_IVch*IBch2
    dDBcl2 = (1-cp)*rhov_IVcl2*gIV_IVcl*IBcl2
    dDBah2 = (1-cp)*rhov_IVah2*gIV_IVah*IBah2
    dDBal2 = (1-cp)*rhov_IVal2*gIV_IVal*IBal2
    dDBeh2 = (1-cp)*rhov_IVeh2*gIV_IVeh*IBeh2
    dDBel2 = (1-cp)*rhov_IVel2*gIV_IVel*IBel2

    res=c(dSch,dScl,dSah,dSal,dSeh,dSel,
          dEch,dEcl,dEah,dEal,dEeh,dEel,
          dAch,dAcl,dAah,dAal,dAeh,dAel,
          dIch,dIcl,dIah,dIal,dIeh,dIel,
          dRch,dRcl,dRah,dRal,dReh,dRel,
          dDch,dDcl,dDah,dDal,dDeh,dDel,
          
          dSch2,dScl2,dSah2,dSal2,dSeh2,dSel2,
          dEch2,dEcl2,dEah2,dEal2,dEeh2,dEel2,
          dAch2,dAcl2,dAah2,dAal2,dAeh2,dAel2,
          dIch2,dIcl2,dIah2,dIal2,dIeh2,dIel2,
          dRch2,dRcl2,dRah2,dRal2,dReh2,dRel2,
          dDch2,dDcl2,dDah2,dDal2,dDeh2,dDel2,
          
          dSVch,dSVcl,dSVah,dSVal,dSVeh,dSVel,
          dEVch,dEVcl,dEVah,dEVal,dEVeh,dEVel,
          dAVch,dAVcl,dAVah,dAVal,dAVeh,dAVel,
          dIVch,dIVcl,dIVah,dIVal,dIVeh,dIVel,
          dRVch,dRVcl,dRVah,dRVal,dRVeh,dRVel,
          dDVch,dDVcl,dDVah,dDVal,dDVeh,dDVel,
          
          dSVch2,dSVcl2,dSVah2,dSVal2,dSVeh2,dSVel2,
          dEVch2,dEVcl2,dEVah2,dEVal2,dEVeh2,dEVel2,
          dAVch2,dAVcl2,dAVah2,dAVal2,dAVeh2,dAVel2,
          dIVch2,dIVcl2,dIVah2,dIVal2,dIVeh2,dIVel2,
          dRVch2,dRVcl2,dRVah2,dRVal2,dRVeh2,dRVel2,
          dDVch2,dDVcl2,dDVah2,dDVal2,dDVeh2,dDVel2,
          
          dSBch,dSBcl,dSBah,dSBal,dSBeh,dSBel,
          dEBch,dEBcl,dEBah,dEBal,dEBeh,dEBel,
          dABch,dABcl,dABah,dABal,dABeh,dABel,
          dIBch,dIBcl,dIBah,dIBal,dIBeh,dIBel,
          dRBch,dRBcl,dRBah,dRBal,dRBeh,dRBel,
          dDBch,dDBcl,dDBah,dDBal,dDBeh,dDBel,

          dSBch2,dSBcl2,dSBah2,dSBal2,dSBeh2,dSBel2,
          dEBch2,dEBcl2,dEBah2,dEBal2,dEBeh2,dEBel2,
          dABch2,dABcl2,dABah2,dABal2,dABeh2,dABel2,
          dIBch2,dIBcl2,dIBah2,dIBal2,dIBeh2,dIBel2,
          dRBch2,dRBcl2,dRBah2,dRBal2,dRBeh2,dRBel2,
          dDBch2,dDBcl2,dDBah2,dDBal2,dDBeh2,dDBel2)
    # cat(sum(res),"\n")
    list(res)
  }) 
} 