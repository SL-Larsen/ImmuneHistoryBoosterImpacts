# Immune history influences SARS-CoV-2 booster impacts: the role of efficacy and redundancy
Code and data for: "Immune history influences SARS-CoV-2 booster impacts: the role of efficacy and redundancy"

All codes and data necessary to produce simulation data, aggregate and analyse this data, and produce the figures for this paper are included. 

## Requirements
This code can be run with R versions 4.2.1-4.4.1 and Python versions 3.8.13-3.12.4.

## Running the Code
Raw simulation data for the History Specific Model is included instead of aggregated simulation data, as these files are extremely large and not storable by GitHub (e.g. individual variant-specific immune histories over time for thousands of simulations). These aggregate data can be generated with "Larsen_et_al_Analysis.Rmd". These data can also be made available to individuals upon reasonable request, but will not be stored in a public repository. 

"antigenicSin" files contain Python code to produce History Specific Model simulation data. These are written to run on a high-performance computer by reading off of an array of simulation parameters which we have included in the "Data" folder. To test this code on your local computer, we recommend that you run only 1 replicate per scenario as the expected run time across all scenarios and replicates included in the paper is several days to 1 week. Folders marked with "ODE SA" contain R code to run the Hybrid Immunity Model, and the expected run time is 20 minutes on your local computer.  Use "Larsen_et_al_Analysis.Rmd" to aggregate simulation data into files for main ("Larsen_et_al_Main.Rmd", "Larsen_et_al_Supplement.Rmd") and supplementary figures. 

## Data Sources
Neutralizing antibody titers data were obtained from [1,2]. SARS-CoV-2 household secondary attack rates where obtained from [3]. Population estimates and projections come from census data [4]. Reported cases were downloaded from the WHO portal [5]. Country-specific contact matrices were obtained from [6]. Estimates of relative mobility by SES where obtained from [7]. Within-SES contact vs. across-SES were inferred using data from [8, 9]. SES-stratified vaccination rates were obtained from [10]. Booster coverage was downloaded from Our World Data website [11]. Monovalent booster efficacy values with respect to the bivalent booster were obtained from [12, 13]. Vaccine efficacy values used to parameterize the hybrid-immunity model were obtained from [14,15, 16, 17]. 

## References

[1] Manali, M., Bissett, L. A., Amat, J. A., Logan, N., Scott, S., Hughes, E. C., ... & Murcia, P. R. (2023). SARS-CoV-2 evolution and patient immunological history shape the breadth and potency of antibody-mediated immunity. The Journal of infectious diseases, 227(1), 40-49.

[2] Suryawanshi, R. K., Chen, I. P., Ma, T., Syed, A. M., Brazer, N., Saldhi, P., ... & Ott, M. (2022). Limited cross-variant immunity from SARS-CoV-2 Omicron without vaccination. Nature, 607(7918), 351-355.

[3] Madewell, Z. J., Yang, Y., Longini, I. M., Halloran, M. E., & Dean, N. E. (2022). Household secondary attack rates of SARS-CoV-2 by variant and vaccination status: an updated systematic review and meta-analysis. JAMA network open, 5(4), e229317-e229317.

[4] "Population estimates and projections for 227 countries and areas." United States Census Bureau, United States Census Bureau, 2022, www.census.gov/data-tools/demo/idb/#/pop?COUNTRY_YEAR=2022&COUNTRY_YR_ANIM=2022&FIPS_SINGLE=EC&FIPS=EC&popPages=BYAGE&POP_YEARS=2022&menu=popViz.

[5] "WHO Coronavirus (COVID-19) Dashboard." , World Health Organization, 2022, covid19.who.int/.

[6] Prem, K., Cook, A. R., & Jit, M. (2017). Projecting social contact matrices in 152 countries using contact surveys and demographic data. PLoS computational biology, 13(9), e1005697.

[7] Mena, G. E., Martinez, P. P., Mahmud, A. S., Marquet, P. A., Buckee, C. O., & Santillana, M. (2021). Socioeconomic status determines COVID-19 incidence and related mortality in Santiago, Chile. Science, 372(6545), eabg5298.

[8] Yechezkel, M., Weiss, A., Rejwan, I., Shahmoon, E., Ben-Gal, S., & Yamin, D. (2021). Human mobility and poverty as key drivers of COVID-19 transmission and control. BMC public health, 21, 1-13.

[9] Bokányi, E., Juhász, S., Karsai, M., & Lengyel, B. (2021). Universal patterns of long-distance commuting and social assortativity in cities. Scientific reports, 11(1), 20829.

[10] Larsen, S. L., Shin, I., Joseph, J., West, H., Anorga, R., Mena, G. E., ... & Martinez, P. P. (2023). Quantifying the impact of SARS-CoV-2 temporal vaccination trends and disparities on disease control. Science Advances, 9(31), eadh9920.

[11] Edouard Mathieu, Hannah Ritchie, Lucas Rodés-Guirao, Cameron Appel, Charlie Giattino, Joe Hasell, Bobbie Macdonald, Saloni Dattani, Diana Beltekian, Esteban Ortiz-Ospina and Max Roser (2020) - "Coronavirus Pandemic (COVID-19)". Published online at OurWorldInData.org. Retrieved from: 'https://ourworldindata.org/coronavirus' [Online Resource]

[12] Lin, D. Y., Xu, Y., Gu, Y., Zeng, D., Wheeler, B., Young, H., ... & Moore, Z. (2023). Effectiveness of bivalent boosters against severe omicron infection. New England Journal of Medicine, 388(8), 764-766.

[13] Mateo-Urdiales, A., Sacco, C., Fotakis, E. A., Del Manso, M., Bella, A., Riccardo, F., ... & Fabiani, M. (2023). Relative effectiveness of monovalent and bivalent mRNA boosters in preventing severe COVID-19 due to omicron BA. 5 infection up to 4 months post-administration in people aged 60 years or older in Italy: a retrospective matched cohort study. The Lancet Infectious Diseases, 23(12), 1349-1359.

[14] Andrews, N., Stowe, J., Kirsebom, F., Toffa, S., Rickeard, T., Gallagher, E., ... & Lopez Bernal, J. (2022). Covid-19 vaccine effectiveness against the Omicron (B. 1.1. 529) variant. New England Journal of Medicine, 386(16), 1532-1546.

[15] "COVID-19 Vaccine AstraZeneca Real-World Evidence Summary." . , Astra Zeneca, 2021, www.astrazeneca.com/content/dam/az/covid-19/media/factsheets/COVID-19_Vaccine_AstraZeneca_Real-World_Evidence_Summary.pdf.

[16] "Boosting with AstraZeneca’s vaccine provides high protection against Omicron, equivalent to mRNA COVID-19 vaccines." Astra Zeneca, Astra Zeneca, 8 Sept. 2022, www.astrazeneca.com/country-sites/thailand/press-releases/boosting-with-astrazenecas-vaccine-provides-high-protection-against-omicron-equivalent-to-mrna-covid-19-vaccines.html.

[17] Solante, R., Alvarez-Moreno, C., Burhan, E., Chariyalertsak, S., Chiu, N. C., Chuenkitmongkol, S., ... & Thwaites, G. (2023). Expert review of global real-world data on COVID-19 vaccine booster effectiveness and safety during the omicron-dominant phase of the pandemic. Expert review of vaccines, 22(1), 1-16.
