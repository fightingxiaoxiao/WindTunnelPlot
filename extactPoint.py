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
            if prop == 'g':
                objectName = line.split()[1]
                if objectName[:6] != 'object':
                    flag = True
                    continue

            if flag:
                points[objectName] = [
                    float(line.split()[1]), float(line.split()[2]), 0.]
                flag = False

    return points


def mergeExcel(filenames):
    """
    合并多个标准测压数据文件.
    """
    dfs = []

    for f in filenames:
        dfs.append(pd.read_excel(f, sheet_name=None, index_col=0))

    fileNum = len(filenames)
    for key in dfs[0]:
        for i in range(fileNum - 1):
            dfs[i + 1][key] = pd.concat([dfs[i][key], dfs[i + 1][key]])
    print(dfs[-1]['0'])

def addToExcel(filename, points, newFilename):
    df_t = pd.read_excel(filename, sheet_name=None, index_col=0)
    writer = pd.ExcelWriter(newFilename)
    for key in df_t:
        df = df_t[key]
        print(key)
        print(df)
        df.insert(0, 'x', 1E32)
        df.insert(1, 'y', 1E32)
        df.insert(2, 'z', 1E32)

        for p_name in points:
            df.loc[p_name, 'x'] = points[p_name][0]
            df.loc[p_name, 'y'] = points[p_name][1]
            df.loc[p_name, 'z'] = points[p_name][2]

        df.to_excel(writer, sheet_name=key)
    writer.close()


def main():
    filename = 'points.obj'
    excelName = 'all243.xlsx'
    newExcelName = 'result.xlsx'
    #mergeExcel(['resultA851.xlsx'])
    points = readPoints(filename)
    addToExcel(excelName, points, newExcelName)


if __name__ == '__main__':
    main()
