#encoding:utf-8
import config
filepath =config.getRecordFundFile()



def recordPathCover(iteration,coverage,uncoverlist):#更新迭代次数、路径覆盖率、未覆盖路径表
    name = "pathCover.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines("第" + str(iteration) + "代 ")
    f.writelines("种群的路径覆盖率: "+str(coverage) + "\n")
    f.writelines("未覆盖的路径： "+str(uncoverlist) + "\n")
    f.write("=====================================\n")
    f.close()


def recordAllM(iteration,M):#记录迭代中所有的评估矩阵信息：每次迭代都不断记录新的评估矩阵
    name = "M.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines("第" + str(iteration) + "代 " + "\n")
    for l in M:
        f.write(str(l) + "\n")
    f.write("=====================================\n")
    f.close()




def recordLocalUserTime(usertime):#记录局部运行总时间
    name = "LocalUserTime.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines("总用时" + str(usertime) + "\n")
    f.write("=====================================\n")
    f.close()


def recordLocalUnCoverPath(ucp):#记录局部未覆盖路径
    name = "UncoverPafterLocal.txt"
    file = filepath + name
    f = open(file, "a+")
    for p in ucp:
        f.writelines( str(p) +"\n")
    f.write("=====================================\n")
    f.close()


def recordTestCase(i,path,test_fit,coverage):#记录测试用例的情况
    name = "RunningTestCase.txt"
    file = filepath + name
    f = open(file, "a+")
    f.write("第"+str(i)+"个测试用例"+str(path) + "\n")
    f.write("个体覆盖路径情况："+str(test_fit) + "\n")
    f.write("个体覆盖率：" + str(coverage) + "\n")
    f.write("================================================" + "\n")
    f.close()


def recordTestCaseRunTime(lenth_ind,time):#记录测试数据运行时间
    name = "TestCaseRunTime.txt"
    file = filepath + name
    f = open(file, "a+")
    f.write(str(lenth_ind)+" "+str(time) + "\n")
    f.close()

def recordPopTime(time):#记录测试数据运行时间
    name = "TestPopRunTime.txt"
    file = filepath + name
    f = open(file, "a+")
    f.write(str(time) + "\n")
    f.close()

def recordCoverPath(pathlist):#记录覆盖路径
    name = "SpathCover.txt"
    file = filepath + name
    f = open(file, "a+")
    f.write(str(pathlist) + "\n")
    f.write("================================================" + "\n")
    f.close()


def recordInitalPop(strf,pop):#记录初始种群
    name = "pop.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines(strf + "\n")
    for op in pop:
        f.write(str(op) + "\n")
    f.write("=====================================\n")
    f.close()

if __name__ == '__main__':
    pathindex = 1
    count = 9
    # recordSAcount(pathindex,count)