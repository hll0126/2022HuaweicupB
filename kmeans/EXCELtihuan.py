import argparse

def main(yuanwenjian, tihuan, result):

    for line1 in open(yuanwenjian):
       # CZXH = 0
        line1 = line1.strip()
        for line2 in open(tihuan):
            target = line2.split(',')[0]
            #print(target)
            if line1 in target:
                #print(target)
                with open(result, 'a') as f:
                    line3 = line2.split(',')[1]
                    #print(line3)
                    #newCZXH = str(CZXH) + '\n'
                    #Wline3 = str(line3) + '\n'
                    f.write(line3)
                    break
            #else:
                #CZXH += 1


main(yuanwenjian='./dingdan22.csv', tihuan='./DDtihuan22.csv', result='.\\DDDpici22.csv')


'''
    with open(result, 'w') as f:


        line_number = 1
        for line in open(tihuan):
            if line_number == 1:
                f.write(line)
                line_number += 1
            else:

                target = line.split('\t')[1]
                if target in Mole_List:
                    f.write(line)
'''