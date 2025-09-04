
# from liste_echantillons.xlsx of Ifremer repository to sample_list.txt 

lsdt = readxl::read_xlsx("liste_echantillons.xlsx")
lsdt = lsdt[-1,] #row of NA
# remove absent samples du to insufficient quality
lsdt = lsdt[!is.element(lsdt$Sample, c("E2267", "E2269", "E2292", "E2296","E2265", "E2270")),] 
write.table(lsdt, file = 'sample_list.txt')

