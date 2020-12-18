import csv
import os
import datetime
import struct
import io

def extractColumns():
	# for convenience, extracting everything as int 
	discountStruct = struct.Struct("i")
	quantityStruct = struct.Struct("i")
	shipDateStruct = struct.Struct("i")
	# output file names
	discountFile = open("tables/discount.col", "wb")
	quantityFile = open("tables/quantity.col", "wb")
	shipDateFile = open("tables/shipdate.col", "wb")
	# specify input files
	with open("tables" + "/" + "lineitem.tbl", newline='\n') as csvfile:
		linereader = csv.reader(csvfile, delimiter='|')
		# skip headers that were manually written
		next(linereader, None)
		for line in linereader:
			# get date
			fulldate = line[10].split("-")
			shipdate = datetime.date(int(fulldate[0]), int(fulldate[1]), \
			int(fulldate[2]))
			# A select query reveals that the minimum date is 1992-01-02
			minShipDate = datetime.date(1992, 1, 2)
			diffShipDate = shipdate - minShipDate
			# when discounts are extracted, there are at most 2 decimal spaces
			discount = int(round(float(line[6]), 2) * 100)
			quantity = int(line[4])
			discountFile.write(discountStruct.pack(discount))
			quantityFile.write(quantityStruct.pack(quantity))
			shipDateFile.write(shipDateStruct.pack(diffShipDate.days))
	discountFile.close()
	quantityFile.close()
	shipDateFile.close()

if __name__ == "__main__":
	extractColumns()