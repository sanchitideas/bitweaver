import csv
import os
import datetime
import struct
import io
if not os.path.exists("tables"):
	os.makedirs("tables")

def extractColumns():
	discountStruct = struct.Struct("B")
	quantityStruct = struct.Struct("B")
	shipDateStruct = struct.Struct("H")
	discountFile = open("tables/arrow_discount.col", "wb")
	quantityFile = open("tables/arrow_quantity.col", "wb")
	shipDateFile = open("tables/arrow_shipdate.col", "wb")
	with open("tables" + "/" + "lineitem.tbl", newline='\n') as csvfile:
		linereader = csv.reader(csvfile, delimiter='|')
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