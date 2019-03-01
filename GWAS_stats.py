from time import time
from scipy.stats import chisquare
import random
start = time()



def main():
    fake_outfile = 'fake_data.tsv'
    true_positives = 1000
    SNPs = 1000000
    sample_size = 500000
    fake_data = create_fake_data(true_positives, SNPs, sample_size)
    fake_FDR_data = FDR(fake_data)
    
    evaluate_FDR(fake_FDR_data, true_positives, SNPs)

    fake_FNR_data = FNR(fake_FDR_data)
    write(fake_FNR_data, fake_outfile)


def GWAS():
    infile = open('Formatted_GWAS.txt', 'r')
    positions = set()
    counter = 0

    for line in infile:
        line = line.split('\t')
        counter += 1
        try:
            positions.add(int(line[1]))
        except:
            pass
    print(counter, len(positions))
    positions = list(positions)
    positions = sorted(positions)
    print(positions[0:10])


def SNPs():
    infile = open('combined_snps.0.3.final.sorted', 'r')
    counter = 0
    genes = set()
    for line in infile:
        line = line.split('\t')
        genes.add(line[4])
        counter += 1
    print(counter)
    print(len(genes))

def create_fake_data(true_positives, SNPs, sample_size):

    SNP_name_SNP_ps = []
    for i in range(SNPs):
        SNP_name = i
        observations = [0,0]
        for j in range(sample_size):
            if i < true_positives:
                x = random.random()
                if x > 0.51:
                    x = True
                else:
                    x = False
            else:
                x = random.choice([True, False])
            if x:
                observations[0] += 1
            else:
                observations[1] += 1
        chi, p = chisquare(observations)
        SNP_name_SNP_ps.append([SNP_name, p])

    return SNP_name_SNP_ps

def takeSecond(elem):
    return elem[1]

def FDR(data):
    data = sorted(data, key=takeSecond, reverse = True)
    FDR_names_values = []
    previous_FDR = 10
    for i in range(len(data)):
        

        if i > 0: 
            predicted_false_discoveries = data[i][1]*len(data)
            FDR = predicted_false_discoveries/(len(data)-(i))
            FDR_names_values.append([data[i][0], data[i][1], min(FDR, previous_FDR)])
            previous_FDR = min(FDR, previous_FDR)
        else:
            FDR_names_values.append([data[i][0], data[i][1], data[i][1]])
            previous_FDR = data[i][1]

    return FDR_names_values





def FNR(data):
    data = sorted(data, key=takeSecond)
    FNR_names_values = []
    previous_FNR = 10
    for i in range(len(data)):
        inverse_p = 1 - data[i][1]
        

        if i > 0: 
            predicted_false_nondiscoveries = inverse_p*len(data)
            FNR = predicted_false_nondiscoveries/(len(data)-(i))
            FNR_names_values.append([data[i][0], data[i][1], data[i][2], min(FNR, previous_FNR)])
            previous_FNR = min(FNR, previous_FNR)
        else:
            FNR_names_values.append([data[i][0], data[i][1], data[i][2], inverse_p])
            previous_FNR = inverse_p

    return FNR_names_values











def evaluate_FDR(data, true_positives, SNPs):
    Rejected_Nulls = 0
    Accepted_Nulls = 0
    Rejected_Positives = 0
    Accepted_Positives = 0
    Total_Null = SNPs-true_positives
    data = sorted(data, key=takeSecond)
    for i in range(len(data)):
        if data[i][0] < true_positives:
            if data[i][2] < 0.05:
                Accepted_Positives += 1
            else:
                Rejected_Positives += 1
        else:
            if data[i][2] < 0.05:
                Accepted_Nulls += 1
            else:
                Rejected_Nulls += 1

    print(Rejected_Nulls, Accepted_Nulls, Rejected_Positives, Accepted_Positives)



def write(data, outfile):
    outfile = open(outfile, 'w')
    outfile.write('{}\t{}\t{}\t{}\n'.format('name', 'p', 'FDR', 'FNR'))
    for i in range(len(data)):
        # print(data[i])
        outfile.write('{}\t{}\t{}\t{}\n'.format(data[i][0], data[i][1], data[i][2], data[i][3]))
        







main()





