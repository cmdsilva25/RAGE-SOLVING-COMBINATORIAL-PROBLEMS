import csv
import sys
from random import randint
import copy
import time

TRACE_FIELD_ID=0
TRACE_FIELD_TIME=1
TRACE_FIELD_CEL=2
TRACE_FIELD_R=3

CANDIDATE_FIELD_ID=0
CANDIDATE_SUM_OF_SCORES=1
CANDIDATE_NUMBER_OF_SCORES=2
CANDIDATE_AVERAGE_SCORE=3
CANDIDATE_ELEMENT_STATUS=4

APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_VEHICLE_ID=0
APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_TOTAL_TRIP=1
APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_DISTANCE=2

STATUS_ELEMENT_NOT_IN_SOLUTION=0
STATUS_ELEMENT_ALREADY_IN_SOLUTION=1
STATUS_ELIMINATED_ELEMENT=2

#OPTIMIZER_TOTAL_NUMBER_OF_CANDIDATE_SOLUTIONS = 10
#OPTIMIZER_INITIAL_NUMBER_OF_FILTERED_SOLUTIONS = 10
#OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION = 2
#OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION = 1
#OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION = 500

OPTIMIZER_TOTAL_NUMBER_OF_CANDIDATE_SOLUTIONS = 10000
#OPTIMIZER_INITIAL_NUMBER_OF_FILTERED_SOLUTIONS = 1000
#OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION = 205
#OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION = 100
#OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION = 1000

def application_dependent_load_data(fname):
    ficheiro = open(fname, 'r')
    reader = csv.reader(ficheiro, delimiter=';')
    trace=[]
    for linha in reader:
        id = int(linha[0])
        time = float(linha[1])
        cel = int(linha[2])*100 + int(linha[3])
        r = int (linha[4])
        item = [int(id),int(time),int(cel),int(r)]
        trace.append(item)
    ficheiro.close()
    return trace

def optimizer_create_candidate_elements(n):
    #each candidate element has the format [id, score_sum, number_of_solutions, average_score, status of the element]
    dict_candidate_elements={}
    for i in range(0,n+1): #sera que eh bug???
        dict_candidate_elements[i]=[i,0,0,0,STATUS_ELEMENT_NOT_IN_SOLUTION]
    #print('dict_candidate_elements:',dict_candidate_elements)
    return dict_candidate_elements

def application_dependent_filter_candidate_elements (dict_candidate_elements, n, trace ):
    for i in trace:
        dict_candidate_elements[ i[TRACE_FIELD_CEL] ][CANDIDATE_AVERAGE_SCORE] += i[TRACE_FIELD_R]

    #sort elements by the score
    list_candidate_elements = sorted(dict_candidate_elements.items(), key=lambda x: x[1][CANDIDATE_AVERAGE_SCORE], reverse=True)

    dict_filtered_candidate_elements={}
    for i in list_candidate_elements[0:n]:
        dict_filtered_candidate_elements[ i[0] ] = i[1]

    #print('dict_filteres_candidate_elements:',len(dict_filtered_candidate_elements),dict_filtered_candidate_elements)
    return dict_filtered_candidate_elements

def optimizer_reduce_set_of_candidate_elements (dict_candidate_elements, n):

    #sort elements by the score
    list_candidate_elements = sorted(dict_candidate_elements.items(), key=lambda x: x[1][CANDIDATE_AVERAGE_SCORE], reverse=True)

    dict_filtered_candidate_elements={}
    for i in list_candidate_elements[0:n]:
        dict_filtered_candidate_elements[ i[0] ] = i[1]

    #print('dict_filteres_candidate_elements:',len(dict_filtered_candidate_elements),dict_filtered_candidate_elements)
    return dict_filtered_candidate_elements


def optimizer_generate_random_solution(dict_candidate_elements, solution_len):
    #each solution has the form: [score]

    keys_list = list( dict_candidate_elements.keys() )

    selected_elements=0
    dict_generated_solution={}
    while selected_elements<solution_len:
        index = randint(0, len(dict_candidate_elements)-1)
        key = int(keys_list[index])
        if dict_candidate_elements[key][CANDIDATE_ELEMENT_STATUS]==STATUS_ELEMENT_NOT_IN_SOLUTION:
            dict_candidate_elements[key][CANDIDATE_ELEMENT_STATUS]=STATUS_ELEMENT_ALREADY_IN_SOLUTION
            selected_elements+=1
            dict_generated_solution[ dict_candidate_elements[key][CANDIDATE_FIELD_ID] ]=0
            #print (selected_elements)

    #reset the status of elements
    for dgs in dict_generated_solution:
        dict_candidate_elements[dgs][CANDIDATE_ELEMENT_STATUS]=STATUS_ELEMENT_NOT_IN_SOLUTION

    #print('generate_random_solution:', len(dict_generated_solution), dict_generated_solution)
    return dict_generated_solution

def application_dependent_compute_score(trace, solution):
    #dict_veiculos has the form [total trip time, trip under coverage]
    dict_veiculos={}
    for item in trace:
        if item[TRACE_FIELD_ID] not in dict_veiculos:
            dict_veiculos[ item[TRACE_FIELD_ID] ] = [item[TRACE_FIELD_ID],0,0]

        if item[TRACE_FIELD_CEL] in solution:
            dict_veiculos[ item[TRACE_FIELD_ID] ] = [
                item[TRACE_FIELD_ID],
                dict_veiculos[ item[TRACE_FIELD_ID] ][APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_TOTAL_TRIP] + item[TRACE_FIELD_R],
                dict_veiculos[ item[TRACE_FIELD_ID] ][APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_DISTANCE] + item[TRACE_FIELD_R]
            ]
        else:
            dict_veiculos[ item[TRACE_FIELD_ID] ] = [
                item[TRACE_FIELD_ID],
                dict_veiculos[ item[TRACE_FIELD_ID] ][APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_TOTAL_TRIP] + item[TRACE_FIELD_R],
                dict_veiculos[ item[TRACE_FIELD_ID] ][APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_DISTANCE]
            ]

    #calcular score
    dict_razao_alpha={}
    for i in range(0,11,1):
        dict_razao_alpha[i]=0

    for k,v in dict_veiculos.items():
        distancia_percorrida = v[APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_TOTAL_TRIP]
        distancia_coberta    = v[APPLICATION_DEPENDENT_DICT_VEHICLES_FIELD_DISTANCE]
        percentual           = int(100*distancia_coberta/distancia_percorrida)

        #print('distancia_percorrida: ', distancia_percorrida)
        #print('distancia_coberta: ', distancia_coberta)

        if percentual>=0:
            dict_razao_alpha[0]+=1
        if percentual>=10:
            dict_razao_alpha[1]+=1
        if percentual>=20:
            dict_razao_alpha[2]+=1
        if percentual>=30:
            dict_razao_alpha[3]+=1
        if percentual>=40:
            dict_razao_alpha[4]+=1
        if percentual>=50:
            dict_razao_alpha[5]+=1
        if percentual>=60:
            dict_razao_alpha[6]+=1
        if percentual>=70:
            dict_razao_alpha[7]+=1
        if percentual>=80:
            dict_razao_alpha[8]+=1
        if percentual>=90:
            dict_razao_alpha[9]+=1
        if percentual>=100:
            dict_razao_alpha[10]+=1

    #calcular score
    score = 0
    for i in range(1,11):
        score += 0.1 * ((dict_razao_alpha[i]+dict_razao_alpha[i-1]) / (2.0*len(dict_veiculos)))

    return score


def application_dependent_save_result(best_solution_ever_seen, best_score_ever_seen, flogname):
    list = best_solution_ever_seen.items()

    # print('creating output file: ', fname)
    fout = open(flogname+'.rsus.csv', "w")
    for l in list:
        fout.write (str(int(l[0]/100)) + ';' + str(l[0]%100) + ';;')
    fout.write('\n' + str(best_score_ever_seen) + '\n')
    fout.close()

def run_simulation(trace, flogname, simulation_count,
                   OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION,
                   OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION,
                   OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION,
                   OPTIMIZER_INITIAL_NUMBER_OF_FILTERED_SOLUTIONS):
    #initial_time = time.process_time()

    #create set of candidate elements to appear in solution number from 0..OPTIMIZER_TOTAL_NUMBER_OF_CANDIDATE_SOLUTIONS-1
    dict_candidate_elements = optimizer_create_candidate_elements(OPTIMIZER_TOTAL_NUMBER_OF_CANDIDATE_SOLUTIONS)

    dict_filtered_candidate_elements = application_dependent_filter_candidate_elements( dict_candidate_elements, OPTIMIZER_INITIAL_NUMBER_OF_FILTERED_SOLUTIONS, trace )

    #start the musical chairs
    number_of_available_candidates = OPTIMIZER_INITIAL_NUMBER_OF_FILTERED_SOLUTIONS

    # reset the score from previous runs
    for dfce in dict_filtered_candidate_elements:
        dict_filtered_candidate_elements[dfce][CANDIDATE_SUM_OF_SCORES] = 0
        dict_filtered_candidate_elements[dfce][CANDIDATE_NUMBER_OF_SCORES] = 0

    #create log file
    fout = open(flogname, "w")
    fout.write ('OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION;' + str(OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION) +'\n')
    fout.write ('OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION;' + str(OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION) +'\n')
    fout.write ('OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION;' + str(OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION) + '\n')
    fout.write ('simulation_cnt;cen;iter;worse;score;best;time\n')

    best_solution_ever_seen={}
    best_score_ever_seen=-1000

    worse_solution_ever_seen={}
    worse_score_ever_seen=10000

    #now we loop reducing the number of candidate elements until reaching the OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION
    while (number_of_available_candidates>OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION):

        #print ('number_of_available_candidates=', number_of_available_candidates)

        #print('iteration:', number_of_available_candidates)
        #in order to remove candidate elements, we must generate several solutions and make a tournament
        #allowing us to infer the individual contribution of each element to the score
        for i in range(0, OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION):
            dict_actual_solution = optimizer_generate_random_solution(dict_filtered_candidate_elements, OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION)
            #print('random solution generated')
            score = application_dependent_compute_score(trace,dict_actual_solution)

            #check if this is the best ever seen score so far
            if score>best_score_ever_seen:
                best_score_ever_seen=score
                best_solution_ever_seen=copy.deepcopy(dict_actual_solution)

            if score<worse_score_ever_seen:
                worse_score_ever_seen=score
                worse_solution_ever_seen=copy.deepcopy(dict_actual_solution)

            if i%1==0:
                #print('actual score:', score, 'best score: ', best_score_ever_seen, 'number of candidates=', number_of_available_candidates)
                #fout.write('%d;%d;%f;%f;%f;%f\n'%(number_of_available_candidates,i,worse_score_ever_seen,score,best_score_ever_seen,time.process_time()-initial_time))
                fout.write('%d;%d;%d;%f;%f;%f\n'%(simulation_count,number_of_available_candidates,i,worse_score_ever_seen,score,best_score_ever_seen))
                fout.flush()

            #update the score in dict_filtered_candidate_elements
            #print ('update the score in dict_filtered_candidate_elements')

            for das in dict_actual_solution:
                dict_filtered_candidate_elements[das][CANDIDATE_SUM_OF_SCORES]+=score
                dict_filtered_candidate_elements[das][CANDIDATE_NUMBER_OF_SCORES]+=1

        #compute average score for each dict_filtered_candidate_elements
        for dfce in dict_filtered_candidate_elements:
            if dict_filtered_candidate_elements[dfce][CANDIDATE_NUMBER_OF_SCORES]!=0:
                dict_filtered_candidate_elements[dfce][CANDIDATE_AVERAGE_SCORE] = dict_filtered_candidate_elements[dfce][CANDIDATE_SUM_OF_SCORES]/dict_filtered_candidate_elements[dfce][CANDIDATE_NUMBER_OF_SCORES]
            else:
                dict_filtered_candidate_elements[dfce][CANDIDATE_AVERAGE_SCORE] = 0

        #sort dict_filtered_candidate_elements according to score
        number_of_available_candidates-=OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION
        dict_filtered_candidate_elements = optimizer_reduce_set_of_candidate_elements(dict_filtered_candidate_elements, number_of_available_candidates)

        print('best score: ', best_score_ever_seen, 'number of candidates=', number_of_available_candidates)

    application_dependent_save_result(best_solution_ever_seen, best_score_ever_seen, flogname)

    #create best score log file for computing stddev
    print("best solution recorded for iteration ", i)
    file_log_best_score = open(flogname + '.best.csv', 'a')
    file_log_best_score.write (str(simulation_count) + ';' + str(best_score_ever_seen) + '\n')
    file_log_best_score.close()
    #return [best_solution_ever_seen, best_score_ever_seen, worse_solution_ever_seen, worse_score_ever_seen]

def main():
    if len(sys.argv)<8:
        print('optimizer <trace> <num_rsus_in_final_solution> <iterations_per_round> <number_of_solutions_eliminated_per_round> <number_of_initial_solutions> <logfile> <repetitions>')
        exit(1)

    trace_file = sys.argv[1]
    OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION = int(sys.argv[2])
    OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION = int(sys.argv[3])
    OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION = int(sys.argv[4])
    OPTIMIZER_INITIAL_NUMBER_OF_FILTERED_SOLUTIONS = int(sys.argv[5])
    flogname = sys.argv[6]
    repetitions = int(sys.argv[7])

    #load data used by optimizer: this is application_dependent
    trace=application_dependent_load_data(trace_file)
    print('trace in memory')

    for simulation_count in range(0,repetitions,1):
        run_simulation(trace, flogname, simulation_count,
        OPTIMIZER_NUMBER_OF_ELEMENTS_IN_FINAL_SOLUTION,
        OPTIMIZER_NUMBER_OF_RANDOM_SOLUTIONS_GENERATED_PER_ITERATION,
        OPTIMIZER_NUMBER_OF_CANDIDATE_SOLUTIONS_ELIMINATED_AFTER_EACH_ITERATION,
        OPTIMIZER_INITIAL_NUMBER_OF_FILTERED_SOLUTIONS)

main()

def imprime_resultado(dict_val):
    list_val = sorted(dict_val.items(), key=lambda x: x[1][3], reverse=True)

    # print('creating output file: ', fname)
    fevolution = open('RESULTADO-OTIMIZACAO.csv', "w")
    fevolution.write('x;y;score;times;st\n')

    # for i in range(0,min(5, len(list_val))):
    for i in range(0, len(list_val)):
        fevolution.write(
            str(int(list_val[i][0] / 100)) + ';' + str(list_val[i][0] % 100) + ';' + str(list_val[i][1][1]) + ';' + str(list_val[i][1][2]) + ';' + str(list_val[i][1][3]) + '\n')
    fevolution.close()

    feliminar = open('RSUS-ELIMINAR.csv', "w")
    for i in range(250, len(list_val)):
        feliminar.write(
            str(int(list_val[i][0])) + '\n')
    feliminar.close()

    fsolucao_atual = open('SOLUCAO-ATUAL-250.csv', "w")
    for i in range(0, 205):
        fsolucao_atual.write(
            str(int(list_val[i][0]/100)) + ';' + str(int(list_val[i][0]%100)) + ';;')
    fsolucao_atual.close()



def especifico_solucao_para_memoria(fname, dps):
    ficheiro = open(fname, 'r')
    reader = csv.reader(ficheiro, delimiter=';')
    dict_solucao={}
    for i in range(0,10000,1):
        #solucao: celula,score,se esta eliminada
        dict_solucao[i]=[i,0,0]

    for linha in reader:
        cont=0
        while cont<dps:
            x = int(linha[cont*3])
            y = int(linha[cont*3+1])
            cel = x*100+y
            dict_solucao[cel]=[cel,0,STATUS_ELEMENT_ALREADY_IN_SOLUTION]
            cont+=1
        ficheiro.close()
        return dict_solucao

