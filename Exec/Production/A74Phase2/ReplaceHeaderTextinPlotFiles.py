import sys

if(len(sys.argv)==1):
 print("Error! No plot file name provided")
 sys.exit()
else:
 with open(sys.argv[1]) as file:
    lines = [line.rstrip() for line in file]


for l in lines:
 filename=l+"/particles/Header"

 search_text = "Version_Two_Dot_One_double"
 replace_text = "Version_Two_Dot_Zero_double"
  
 with open(filename, 'r') as file:
  data = file.read()
  data = data.replace(search_text, replace_text)
  
 with open(filename, 'w') as file:
  file.write(data)
  
