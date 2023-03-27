import strutils

type Record = object
  id, ix: int
  data: seq[float]

let f = open("clust-90-indel-16-samples-152-include-5_K_9_R_1")

while true:
  let line = f.readLine
  if line == "Admixture analysis":
    discard f.readLine()
    discard f.readLine()
    discard f.readLine()
    break

var records: seq[Record]
while true:
  let line = f.readLine
  if line.isEmptyOrWhitespace:
    break
  let 
    colonSplit = line.split(":")
    meta = colonSplit[0].splitWhitespace
    dataStrings = colonSplit[1].splitWhitespace
    id = meta[2]
    ix = meta[4] 
  var data = newSeq[float](dataStrings.len) 
  for i in 0 ..< data.len:
    data[i] = parseFloat(dataStrings[i])
  records.add(Record(id:id, ix:ix, data:data))

f.close()