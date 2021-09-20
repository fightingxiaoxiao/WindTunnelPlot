import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd
import os


def convertRawDataFromExcelToCSV(filename, angle, valueType, region_head):
    data = pd.read_excel(filename, sheet_name=str(angle))
    region_name = []
    count = 0
    for n in data['name']:
        if n.split('-')[0] == region_head:
            region_name.append(count)
        count += 1
    data = data.loc[:, (['name', 'x', 'y', 'z', valueType])]
    data = data.loc[region_name]
    # print(data)
    if not os.path.exists('./contour_data'):
        os.mkdir('./contour_data')

    data.to_csv('./contour_data/'+region_head +
                '_'+str(angle)+'_'+valueType+'.csv', header=0, index=0)


def readCSV(filename):
    """
    读取数据点.
    """
    x, y, z = [], [], []
    with open(filename, 'r') as f:
        for line in f:
            dat = line.split(',')
            x.append(float(dat[0]))
            y.append(float(dat[1]))
            z.append(float(dat[3]))
    return x, y, z
# First create the x and y coordinates of the points.


def readCellInfo(filename):
    tri = []
    with open(filename, 'r') as f:
        for line in f:
            sub_tri = []
            nodes = line.split(',')
            for n in nodes:
                sub_tri.append(int(n))
            tri.append(sub_tri)
    return np.array(tri)


def contour(dataFile, cellFile, region_name, angle, value_type, delta_level, saveFileFolder=None, coeff=1):
    x, y, z = readCSV(dataFile)
    # triang = tri.Triangulation(x, y)
    triang = readCellInfo(cellFile)
    fig1, ax1 = plt.subplots()

    z = np.array(z)
    z *= coeff
    #print(delta_level)
    up = int(np.max(z)/delta_level)
    down = int(np.min(z)/delta_level)

    if up > 0:
        up += 1

    if down < 0:
        down -= 1

    level = [delta_level*i for i in range(down, up)]
    ax1.set_aspect('equal')
    tcf = ax1.tricontourf(x, y, triang, z, cmap='jet', levels=level)
    fig1.colorbar(tcf)
    ax1.tricontour(x, y, triang, z, colors='k',
                   linestyles='solid', linewidths=1, levels=level)
    ax1.set_title(region_name+'_'+str(angle)+'_'+value_type)
    plt.axis('off')

    if saveFileFolder is None:
        plt.show()
    else:
        saveFileName = region_name+'_'+str(angle)+'_'+value_type+'.jpg'
        if not os.path.exists('./'+saveFileFolder):
            os.mkdir('./'+saveFileFolder)
        plt.savefig('./'+saveFileFolder+'/'+saveFileName, dpi=400)


def processList(filename, region_list, angle_list, type_list, iters, work_dir, delta_level, saveFileFolder, coeff=1):
    for r in region_list:
        for a in angle_list:
            for t in type_list:
                processSingle(filename, r, a, t, iters,
                              work_dir, delta_level, saveFileFolder, coeff)


def processSingle(filename, region, angle, value_type, iters, work_dir, delta_level, saveFileFolder=None, coeff=1, platform='windows',block_x=100, block_y=100):
    os.chdir(work_dir)
    dataFile = work_dir+'/out_tmp/'+region + \
        '_'+str(angle)+'_'+value_type+'_interp'+'.csv'
    cellFile = work_dir+'/out_tmp/'+region+'_' + \
        str(angle)+'_'+value_type+'_cell_info'+'.csv'
    zoneAvgFile = work_dir+'/out_coeff/'+region+'_' + \
        str(angle)+'_'+value_type+'_coeff'+'.csv'

    if not os.path.exists('./out_tmp'):
        os.mkdir('./out_tmp')
    if not os.path.exists('./out_coeff'):
        os.mkdir('./out_coeff')

    convertRawDataFromExcelToCSV(filename, angle, value_type, region)
    if platform=="windows":
        execute = 'wsl -e ./windtunnel_tool'
    elif platform=="linux":
        execute = './windtunnel_tool'
    command = execute+' ' + \
        work_dir+'/data/'+region+'.obj' + ' ' + \
        work_dir+'/contour_data/'+region+'_'+str(angle)+'_'+value_type+'.csv' + ' ' + \
        dataFile + ' ' + \
        cellFile + ' ' + \
        zoneAvgFile + ' ' + \
        str(iters) + ' ' + \
        str(block_x) + ' ' + \
        str(block_y)
    # print(command)
    os.system(command)
    contour(dataFile, cellFile, region,
            angle, value_type, delta_level, saveFileFolder, coeff)


if __name__ == "__main__":
    """
    value_type = 'avg'  # avg: 平均值, max: 最大值, min: 最小值, std: 标准差
    region_name = 'AC'
    angle = 270
    iteration = 100000    #work_dir = 'F:/wind_tunnel_data'
    work_dir = '.'
    filename = 'out_2.xlsx'

    processSingle(filename, region_name, angle, value_type,
                  iteration, work_dir, 0.1, saveFileFolder=None, coeff=1/0.5/1.225/6/6, platform="linux",block_x=300, block_y=100)
    """
    
    type_list = ['avg']
    region_list = ['AU', 'BU', 'CU']
    angle_list = [i*15 for i in range(24)]
    iteration = 5000
    work_dir = 'F:/wind_tunnel_data'
    filename = 'out.xlsx'
    processList(filename, region_list, angle_list, type_list,
                iteration, work_dir, 0.2, 'contour', 1/0.5/1.225/6/6, platform="linux")

