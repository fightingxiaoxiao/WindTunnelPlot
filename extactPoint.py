from os import read
import pandas as pd


def readPoints(filename):
    points = {}
    objectName = None
    with open(filename, 'r') as f:
        flag = False
        for line in f:
            line = line[:-1]
            try:
                prop = line.split()[0]
            except IndexError:
                pass
            if prop == 'o':
                objectName = line.split()[1]
                if objectName[:6] != 'object':
                    flag = True
                    continue

            if flag:
                points[objectName] = [
                    float(line.split()[1]), float(line.split()[2]), 0.]
                flag = False

    return points


def addToExcel(filename, points):
    df_t = pd.read_excel(filename, sheet_name=None, index_col=0)
    for key in df_t:
        df = df_t[key]
        df.insert(0, 'x', 10000)
        df.insert(1, 'y', 10000)
        df.insert(2, 'z', 10000)

    for p_name in points:
        df.loc[p_name,'x'] = points[p_name][0]
        df.loc[p_name,'y'] = points[p_name][1]
        df.loc[p_name,'z'] = points[p_name][2]
        #print(df.loc[p_name])
    
    df.to_excel('test.xlsx',sheet_name='0',index='name')

def main():
    filename = 'TEST01.obj'
    excelName = 'resultA851.xlsx'
    points = readPoints(filename)
    addToExcel(excelName, points)


if __name__ == '__main__':
    main()
