#encoding:utf-8
from __future__ import division
import sys
sys.path.append("..")
sys.path.append("../..")
sys.path.append("../module")
from datetime import datetime
import obtain_efsm_info2    # 模型处理器
import generate_seq_ga_check  # 序列生成器-----产生覆盖敏感路径的序列
import execute
import local_search_alg
import recordFun
import config
def main():
    fileName_F = 'D:\programs\PyCharm 2020.1.2\projects\graphTraversal\dataset'
    for L in ["\M.txt","\pathCover.txt","\pop.txt","\RunningTestCase.txt","\SpathCover.txt","\StatisticResult.txt","\TestCaseRunTime.txt","\TestPopRunTime.txt"]:
        fileName = fileName_F + L;
        f = open(fileName, "r+")
        f.truncate()
        f.close()
    SM = obtain_efsm_info2.obtain_efsm()
    popsize,pc,pm,Max = config.getPopParameter()
    for i in range(0, 1):
        print "新一轮开始：", i
        FlagLunhui = True
        startTime = datetime.now()
        pop_time0=datetime.now()
        #pop = generate_seq_ga_check.initialpop_feasible(popsize, SM)
        pop = [['T2', 'T4'],
        ['T1', 'T15', 'T10'],
        ['T2', 'T4', 'T3'],
        ['T1', 'T11', 'T28']]


        strf = "初始种群"
        recordFun.recordInitalPop(strf,pop)#构造初始种群
        print ("种群为： ")
        for sub in pop:
            print (sub)
        #recordFun.recordTestCaseRunTime(i)
        recordFun.recordTestCase(i, 0, 0, 0)
        testdata, M,V = execute.gdata(pop, SM)  # 为序列产生数据，执行测试用例，获取评估矩阵M
        SVTCM, ICM = generate_seq_ga_check.SVTCM_ICM(pop, M, V, 1,[],[])  # 测试用例对耦合矩阵、个体被耦合矩阵
        print ("通过M获得所有敏感路径入口点覆盖标识cover_flag及敏感路径完全覆盖标识flag")
        flag, cover_flag = execute.table_handle(M)
        print ("1敏感路径入口点覆盖标识cover_flag = ",cover_flag)
        print ("1敏感路径完全覆盖标识flag = ", flag)
        coverage, Uncover_path_list = execute.table_coverage(M)#M是评估矩阵，基于M获取违背覆盖的敏感路径集Uncover——path——list）
        ######################## 以上为评估初始种群############################
        iteration = 0#迭代次数
        pop_time1 = datetime.now()
        pop_time = pop_time1 - pop_time0
        recordFun.recordPopTime(pop_time)
        while True:
            recordFun.recordPathCover(iteration, coverage, Uncover_path_list)#路径更新未覆盖路径表
            recordFun.recordAllM(iteration,M)#更新评估矩阵
            pop_time0 = datetime.now()
            if iteration > Max:
                break
            # 先判断种群是否覆盖所有入口点
            print ("敏感路径入口点覆盖标识cover_flag = ", cover_flag)
            print ("敏感路径完全覆盖标识flag = ", flag)
            if cover_flag == 1:  # 如果覆盖了所有敏感路径的入口点
                if flag == 1:  # 在判断是否覆盖所有的敏感路径
                    testcase = zip(pop, testdata)  # 测试用例等于序列+数据
                    print (testcase)
                    break
                else:  # 进行局部搜索
                    print ("第",iteration , "代" + "进行局部搜索")
                    print ("进入局部搜索的pop", pop)
                    print ("进入局部搜索的testdata", testdata)
                    print ("进入局部搜索的M", M)
                    print ("进入局部搜索的Uncover_path_list", Uncover_path_list)
                    # 局部算法接口在这调用
                    #overgame,testdata = local_search_alg.Local_search_for_SA(pop,testdata,M,Uncover_path_list,SM)
                    overgame,testdata = local_search_alg.Local_search_for_HC(pop,testdata,M,Uncover_path_list,SM)
                    #overgame, testdata = local_search_alg.Local_search_for_GA(pop, testdata, M, Uncover_path_list, SM)#SA局部搜索，跟新Uncover_path_list.
                    if overgame == 1:
                        testcase = zip(pop, testdata)  # 测试用例等于序列+数据
                        print (testcase)
                        break # 退出while ，执行下面的时间记录代码
                    else:
                        tip = "搜索失败，结束"
                        print (tip)
                        return tip  #退出整个程序，不会执行下面的时间记录代码了
            else:
                print ("继续全局GA，搜序列")
                #pop,V,M,testdata,cflag,SVTCM, ICM = generate_seq_ga_check.GA(pop, pc, pm, popsize, SM, M,testdata,V,SVTCM, ICM,iteration)
                pop,V,M,testdata,cflag= generate_seq_ga_check.GA(iteration,pop, pc, pm, popsize, SM, M, testdata, V)  # 原始串行GA
                FlagLunhui = False
                if cflag == 1:
                    testcase = zip(pop, testdata)  # 测试用例等于序列+数据
                    print (testcase)
                    break
                else:
                    flag, cover_flag = execute.table_handle(M)
                    coverage, Uncover_path_list  = execute.table_coverage(M)
                iteration += 1
            pop_time1 = datetime.now()
            pop_time = pop_time1 - pop_time0
            recordFun.recordPopTime(pop_time)
            #得到新的flag, cover_flag后，这里会直接到while True那句话再开始
        endTime = datetime.now()
        usertime = endTime - startTime
        print usertime
        recordFun.recordPopTime(usertime)



if __name__ == '__main__':  # not execute when import as a module
    main()
    # table = []
    # print (len(table))
    # t = [ 2, 0.0, 1.004950495049505, 0.0, 0.8544139861204133]
    # for i in range(2):
    #     table.append(t)
    # print (table)
    """
    ['T1', 'T5', 'T6', 'T7'],
    ['T1', 'T5', 'T27', 'T26'],
    ['T1', 'T12', 'T30', 'T26'],
    ['T2', 'T4', 'T3', 'T8'],
    ['T1', 'T11', 'T28', 'T27', 'T28'],
    ['T1', 'T13', 'T35', 'T36', 'T35'],
    ['T1', 'T12', 'T18', 'T21', 'T12'],
    ['T1', 'T5', 'T27', 'T28', 'T6', 'T7'],
    ['T1', 'T12', 'T30', 'T28', 'T27', 'T26'],
    ['T1', 'T9', 'T33', 'T20', 'T14', 'T11', 'T26'],
    ['T1', 'T5', 'T27', 'T29', 'T23', 'T24', 'T10'],
    ['T1', 'T11', 'T29', 'T31', 'T32', 'T30', 'T29'],
    ['T1', 'T15', 'T12', 'T18', 'T21', 'T13', 'T19', 'T9'],
    ['T1', 'T11', 'T28', 'T27', 'T29', 'T18', 'T21', 'T13']]
    ['T1', 'T48', 'T30', 'T28'],"""
