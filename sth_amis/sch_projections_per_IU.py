from sch_simulation.helsim_RUN_KK import *
from sch_simulation.helsim_FUNC_KK import *
from sch_simulation.helsim_FUNC_KK.results_processing import getBurdens
import pickle
import numpy as np
import time
import os
from multiprocessing import Pool

# should be haematobium, mansoni_low_burden or mansoni_high_burden
species = "haematobium"

IU_SLURM = os.getenv("SLURM_ARRAY_TASK_ID")
#IU_SLURM = 2
num_cores = 10


start = time.time()

'''
    extra functions to output the prevalence data
'''
def adjustR0AndKParams(parameters, R0, k):
    parameters.R0 = R0
    parameters.k = k

    # configure the parameters
    parameters = configure(parameters)
    parameters.psi = getPsi(parameters)
    parameters.equiData = getEquilibrium(parameters)
    return parameters

def generateSimData(seed, params, R0, k):
    np.random.seed(seed)
    # adjust parameters here so that we will use the correct set of parameters 
    # to generate the new population for each simulation
    params = adjustR0AndKParams(params, R0, k)
    # return the starting simData
    return setupSD(params)
    
def constructNTDMCResults(params, res, startYear):
    """
    This function will return a data frame forNTDMC containing appropriate data

    Parameters
    ----------
    params: Parameters
        dataclass containing the parameter names and values;
    res: 
        the results from the simulation along with the state of the population
    Returns
    -------
     Data frame for NTDMC which contains simpler
    data, to do with prevalence of SAC and the whole population.
    """

    # unpack results from the simulation here. These will be called to extract data from later
    results = [ item[ 0 ] for item in res ]
    allSD = [ item[ 1 ] for item in res ]
    for i in range(len(results)):

        singleSimResult = results[i]
        SD = allSD[i]

        _, _, _, dfSAC = returnNTDMCOutputs(params, singleSimResult, [5, 15], "SAC", startYear)
        _, _, _, dfAll = returnNTDMCOutputs(params, singleSimResult, [0,100], "Whole Population", startYear)

        if i == 0:
            NTDMC = pd.concat([dfSAC, dfAll], ignore_index = True)
            NTDMC["draw_0"] = NTDMC['draw_1'].values
            NTDMC = NTDMC.drop('draw_1', axis=1)
        else:
            colname = "draw_" + str(i)

            newColsNTDMC = pd.concat([dfSAC, dfAll], ignore_index = True)
            NTDMC[colname] = newColsNTDMC['draw_1'].values

    return NTDMC


def returnNTDMCOutputs(params, results, ageBand, PopType, startYear, prevThreshold = 0.02, surveyType = "KK2", numReps = 1, nSamples = 2, sampleSize = 100):
    output = extractHostData([results])
    prevalence, _, medium_prevalence, heavy_prevalence, _ = getBurdens(output, params, 
                                                                       numReps, ageBand,
                                                                       params.Unfertilized, 
                                                                       surveyType, nSamples, 
                                                                       sampleSize
                                                                       ) 
    allTimes = output[0].timePoints + startYear
    if PopType == "SAC":
        prevBelowThreshold = (medium_prevalence + heavy_prevalence) < prevThreshold
    newrows = pd.DataFrame(
                        {
                            "year_id": allTimes,
                            "age_start": np.repeat(ageBand[0], len(allTimes)),
                            "age_end": np.repeat(ageBand[1], len(allTimes)),
                            "intensity": np.repeat("None", len(allTimes)),
                            "species": np.repeat(params.species, len(allTimes)),
                            "measure": np.repeat("Prevalence " + PopType, len(allTimes)),
                            "draw_1": prevalence,
                        }
                    )
    df1 = newrows
    newrows = pd.DataFrame(
                        {
                            "year_id": allTimes,
                            "age_start": np.repeat(ageBand[0], len(allTimes)),
                            "age_end": np.repeat(ageBand[1], len(allTimes)),
                            "intensity": np.repeat("None", len(allTimes)),
                            "species": np.repeat(params.species, len(allTimes)),
                            "measure": np.repeat("Medium + Heavy Prevalence " + PopType, len(allTimes)),
                            "draw_1": medium_prevalence + heavy_prevalence,
                        }
                    )
    df1 = pd.concat([df1, newrows], ignore_index = True)
    if PopType == "SAC":
        newrows = pd.DataFrame(
                            {
                                "year_id": allTimes,
                                "age_start": np.repeat(ageBand[0], len(allTimes)),
                                "age_end": np.repeat(ageBand[1], len(allTimes)),
                                "intensity": np.repeat("None", len(allTimes)),
                                "species": np.repeat(params.species, len(allTimes)),
                                "measure": np.repeat("Below  EPHP threshold", len(allTimes)),
                                "draw_1": prevBelowThreshold,
                            }
                        )
        df1 = pd.concat([df1, newrows], ignore_index = True)
    return prevalence, medium_prevalence, heavy_prevalence, df1


############################################################################################################
############################################################################################################
############################################################################################################

'''
    Define names of files and other information for the runs. 
    This is the section which needs the most editing for each IU
'''

if species == "haematobium":
    species_prefix = "Haema_"
    pathCountry = '../Maps-SCH/table_iu_idx_haematobium.csv'
if species == "mansoni_high_burden":
    species_prefix = "Man_High_"
    pathCountry = '../Maps-SCH/table_iu_idx_mansoni.csv'
if species == "mansoni_low_burden":
    species_prefix = "Man_Low_"    
    pathCountry = '../Maps-SCH/table_iu_idx_mansoni.csv'

df_IU_country = pd.read_csv(pathCountry)
iu = df_IU_country['IU_CODE'].values[int(IU_SLURM)]
country = df_IU_country['country'].values[int(IU_SLURM)]

print(country)
print(species)
print(str(iu).zfill(5))
print(species_prefix)


''' 
    change to be whatever the scenario we are running is
'''
if species == "haematobium":
    print("Making projections for haematobium for IU " + str(iu).zfill(5) + "...")
    coverageFileName = 'endgame_inputs_haema/InputMDA_MTP_projections_' + str(iu) + '.xlsx'
else:
    print("Making projections for " + str(species) + " for IU " + str(iu).zfill(5) + "...")
    coverageFileName = 'endgame_inputs_mansoni/InputMDA_MTP_projections_' + str(iu) + '.xlsx'

'''
    Change to whatever demography was used to run the fitting
'''
demogName = 'UgandaRural'

# file name to store coverage information in 
coverageTextFileStorageName = 'Man_MDA_vacc/Man_MDA_vacc_' + str(species) + '_' + str(iu) + '.txt'

'''
    Change to whatever parameter file you are using for these runs
'''
paramFileName = 'SCH_params/' + str(species) + '_params_projections.txt'       

''' 
    this should be the time point that the simulations should begin from.
'''
startYear = 1985

'''
    file name for IU specific parameters
'''
RkFilePath = '../post_AMIS_analysis/InputPars_MTP_' + str(species) + '/InputPars_MTP_' + str(iu) + '.csv'

'''
    numSims should be set to 200
'''
numSims = 200

############################################################

'''
    Read in files and prepare the necessary data for the runs
'''

# read in parameter and coverage files
_ = parse_coverage_input(coverageFileName,
     coverageTextFileStorageName)
# initialize the parameters
params = loadParameters(paramFileName, demogName)
# add coverage data to parameters file
params = readCoverageFile(coverageTextFileStorageName, params)
# add vector control data to parameters
params = parse_vector_control_input(coverageFileName, params)

# read in fitted parameters for IU
simparams = pd.read_csv(RkFilePath)
simparams.columns = [s.replace(' ', '') for s in simparams.columns]

# define the lists of random seeds, R0 and k
seed = simparams.iloc[:, 1].tolist()
R0 = simparams.iloc[:, 2].tolist()
k1 = simparams.iloc[:, 3].tolist()

# setup the output times to be every year
params.outTimings =  np.arange(0, params.maxTime, 1) 

# set the survey type to Kato Katz with duplicate slide
surveyType = "KK2"

paramsToAlter = copy.deepcopy(params)

############################################################

'''
    run the simulations
'''
res = Parallel(n_jobs=num_cores)(
             delayed(doRealizationSurveyCoveragePickle)(adjustR0AndKParams(paramsToAlter, R0[i], k1[i]), surveyType, generateSimData(seed[i], paramsToAlter, R0[i], k1[i])) for i in range(numSims)
         )

simData=[item[1] for item in res]

'''
    We want to make another set of pickle files don't we?
    If so here we specify the name of the pickle file 
'''

# want outputs like <ascaris-folder>/AGO/AGO02049/Asc_AGO02049.p
newOutputSimDataFilePath = f'projections/{species}/{country}/{country}' + str(iu).zfill(5) + f'/{species_prefix}{country}' + str(iu).zfill(5) + '.p'
print("Pickle file name:")
print(newOutputSimDataFilePath)
pickle.dump( simData, open( newOutputSimDataFilePath, 'wb' ) )

NTDMC = constructNTDMCResults(params, res, startYear)    
PrevDatasetFilePath = f'projections/{species}/{country}/{country}' + str(iu).zfill(5) + f'/PrevDataset_{species_prefix}{country}' + str(iu).zfill(5) + '.csv'
print("PrevDataset_species_iu.csv file name:")
print(PrevDatasetFilePath)
NTDMC.to_csv(PrevDatasetFilePath, index=False)

print('Finished projections for ' + str(species) + ' in IU ' + str(iu) + ".")

#################################################################################################################

end = time.time()
elapsed_time = end-start
print("elapsed time in seconds: " + str(elapsed_time) )

