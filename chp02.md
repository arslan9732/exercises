## Exercises and solutions for Chapter 2

### Computations in R

1. Sum 2 and 3 using the `+` operator. [Difficulty: **Beginner**]

**solution:**
  ```{r}
2+3
```

2. Take the square root of 36, use `sqrt()`. [Difficulty: **Beginner**]

**solution:**
  ```{r}
sqrt(36)
```

3. Take the log10 of 1000, use function `log10()`. [Difficulty: **Beginner**]

**solution:**
  ```{r}
log10(1000)
```

4. Take the log2 of 32, use function `log2()`. [Difficulty: **Beginner**]

**solution:**
  ```{r}
log2(32)
```

5. Assign the sum of 2,3 and 4 to variable x. [Difficulty: **Beginner**]

**solution:**
  ```{r}
x =  2+3+4
x <- 2+3+4
```

6. Find the absolute value of  the expression `5 - 145` using the `abs()` function. [Difficulty: **Beginner**]

**solution:**
  ```{r}
abs(5-145)
```

7. Calculate the square root of 625, divide it by 5, and assign it to variable `x`.Ex: `y= log10(1000)/5`, the previous statement takes log10 of 1000, divides it by 5, and assigns the value to variable y. [Difficulty: **Beginner**]

**solution:**
  ```{r}
x = sqrt(625)/5

```

8. Multiply the value you get from previous exercise by 10000, assign it to variable x
Ex: `y=y*5`, multiplies `y` by 5 and assigns the value to `y`. 
**KEY CONCEPT:** results of computations or arbitrary values can be stored in variables we can re-use those variables later on and over-write them with new values. 
[Difficulty: **Beginner**]

**solution:**
  ```{r}
x2 = x*10000
```

### Data structures in R 


10. Make a vector of 1,2,3,5 and 10 using `c()`, and assign it to the `vec` variable. Ex: `vec1=c(1,3,4)` makes a vector out of 1,3,4. [Difficulty: **Beginner**]

**solution:**
  ```{r}
c(1:5,10)
vec1=c(1,2,3,4,5,10)
```

11. Check the length of your vector with length(). 
Ex: `length(vec1)` should return 3. [Difficulty: **Beginner**]

**solution:**
  ```{r}
length(vec1)
```

12. Make a vector of all numbers between 2 and 15.
Ex: `vec=1:6` makes a vector of numbers between 1 and 6, and assigns it to the `vec` variable. [Difficulty: **Beginner**]

**solution:**
  ```{r}
vec=2:15 

```

13. Make a vector of  4s repeated 10 times using the `rep()` function. Ex: `rep(x=2,times=5)` makes a vector of 2s repeated 5 times. [Difficulty: **Beginner**]

**solution:**
  ```{r}
rep(x=4,times=10)
rep(4,10) 
```

14. Make a logical vector with TRUE, FALSE values of length 4, use `c()`.
Ex: `c(TRUE,FALSE)`. [Difficulty: **Beginner**]

**solution:**
  ```{r}
c(TRUE,FALSE,FALSE,TRUE,FALSE)
c(TRUE,TRUE,FALSE,TRUE,FALSE)

```

15. Make a character vector of the gene names PAX6,ZIC2,OCT4 and SOX2.
Ex: `avec=c("a","b","c")` makes a character vector of a,b and c. [Difficulty: **Beginner**]

**solution:**
  ```{r}
c("PAX6","ZIC2","OCT4","SOX2")

```

16. Subset the vector using `[]` notation, and get the 5th and 6th elements. 
Ex: `vec1[1]` gets the first element. `vec1[c(1,3)]` gets the 1st and 3rd elements. [Difficulty: **Beginner**]

**solution:**
  ```{r}
vec1[c(5,6)]

```

17. You can also subset any vector using a logical vector in `[]`. Run the following:
  
  ```{r subsetLogicExercise, eval=FALSE}
myvec=1:5
# the length of the logical vector 
# should be equal to length(myvec) 
myvec[c(TRUE,TRUE,FALSE,FALSE,FALSE)] 
myvec[c(TRUE,FALSE,FALSE,FALSE,TRUE)]
``` 
[Difficulty: **Beginner**]


**solution:**
  ```{r}
myvec=1:5
myvec[c(TRUE,TRUE,FALSE,FALSE,FALSE)] # the length of the logical vector should be equal to length(myvec) 
myvec[c(TRUE,FALSE,FALSE,FALSE,TRUE)] 
```

18. `==,>,<, >=, <=` operators create logical vectors. See the results of the following operations:
  
  ```{r,eval=FALSE}
myvec > 3
myvec == 4
myvec <= 2
myvec != 4
```
[Difficulty: **Beginner**]

19.  Use the `>` operator in `myvec[ ]` to get elements larger than 2 in `myvec` which is described above. [Difficulty: **Beginner**]


**solution:**
  ```{r}
myvec[ myvec > 2 ] 

```


20. Make a 5x3 matrix (5 rows, 3 columns) using `matrix()`.
Ex: `matrix(1:6,nrow=3,ncol=2)` makes a 3x2 matrix using numbers between 1 and 6. [Difficulty: **Beginner**]

**solution:**
  ```{r}
mat=matrix(1:15,nrow=5,ncol=3)

```

21. What happens when you use `byrow = TRUE` in your matrix() as an additional argument?
  Ex: `mat=matrix(1:6,nrow=3,ncol=2,byrow = TRUE)`. [Difficulty: **Beginner**]

**solution:**
  ```{r}
mat=matrix(1:15,nrow=5,ncol=3,byrow = TRUE)

```

22. Extract the first 3 columns and first 3 rows of your matrix using `[]` notation. [Difficulty: **Beginner**]


**solution:**
  ```{r}
mat[1:3,1:3]

```

23. Extract the last two rows of the matrix you created earlier.
Ex: `mat[2:3,]` or `mat[c(2,3),]` extracts the 2nd and 3rd rows.
[Difficulty: **Beginner**]

**solution:**
  ```{r}
mat[4:5,] 
mat[c(nrow(mat)-1,nrow(mat)),] 
tail(mat,n=1)
tail(mat,n=2) 
```

24. Extract the first two columns and run `class()` on the result.
[Difficulty: **Beginner**]

**solution:**
  ```{r}
class(mat[,1:2])

```

25. Extract the first column and run `class()` on the result, compare with the above exercise.
[Difficulty: **Beginner**]

**solution:**
  ```{r}
class(mat[,1])

```

26. Make a data frame with 3 columns and 5 rows. Make sure first column is a sequence
of numbers 1:5, and second column is a character vector.
Ex: `df=data.frame(col1=1:3,col2=c("a","b","c"),col3=3:1) # 3x3 data frame`.
Remember you need to make a 3x5 data frame. [Difficulty: **Beginner**]

**solution:**
  ```{r}
df=data.frame(col1=1:5,col2=c("a","2","3","b","c"),col3=5:1)

```

27. Extract the first two columns and first two rows. 
**HINT:** Use the same notation as matrices. [Difficulty: **Beginner**]

**solution:**
  ```{r}
df[,1:2]

df[1:2,]
df[1:2,1:2] 
```

28. Extract the last two rows of the data frame you made.
**HINT:** Same notation as matrices. [Difficulty: **Beginner**]

**solution:**
  ```{r}
df[,4:5]

```

29. Extract the last two columns using the column names of the data frame you made. [Difficulty: **Beginner**]


**solution:**
  ```{r}
df[,c("col2","col3")]

```

30. Extract the second column using the column names.
You can use `[]` or `$` as in lists; use both in two different answers. [Difficulty: **Beginner**]

**solution:**
  ```{r}
df$col2
df[,"col2"]
class(df["col1"])
class(df[,"col1"]) 
```

31. Extract rows where the 1st column is larger than 3.
**HINT:** You can get a logical vector using the `>` operator
, and logical vectors can be used in `[]` when subsetting. [Difficulty: **Beginner**]

**solution:**
  ```{r}
df[df$col1 >3 , ]

```

32. Extract rows where the 1st column is larger than or equal to 3.
[Difficulty: **Beginner**]

**solution:**
  ```{r}
df[df$col1 >= 3 , ]

```

33. Convert a data frame to the matrix. **HINT:** Use `as.matrix()`.
Observe what happens to numeric values in the data frame. [Difficulty: **Beginner**]

**solution:**
  ```{r}
class(df[,c(1,3)])
as.matrix(df[,c(1,3)])
as.matrix(df) 
```


34. Make a list using the `list()` function. Your list should have 4 elements;
the one below has 2. Ex: `mylist= list(a=c(1,2,3),b=c("apple,"orange"))`
[Difficulty: **Beginner**]

**solution:**
  ```{r}
mylist= list(a=c(1,2,3),
             b=c("apple","orange"),
             c=matrix(1:4,nrow=2),
             d="hello") 
```

35. Select the 1st element of the list you made using `$` notation.
Ex: `mylist$a` selects first element named "a".
[Difficulty: **Beginner**]

**solution:**
  ```{r}
mylist$a

```

36. Select the 4th element of the list you made earlier using `$` notation. [Difficulty: **Beginner**]

**solution:**
  ```{r}
mylist$d

```


37. Select the 1st element of your list using `[ ]` notation. 
Ex: `mylist[1]` selects the first element named "a", and you get a list with one element. `mylist["a"]` selects the first element named "a", and you get a list with one element.
[Difficulty: **Beginner**]

**solution:**
  ```{r}
mylist[1] # -> still a list
mylist[[1]] # not a list

mylist["a"] 
mylist[["a"]] 

```
38. Select the 4th element of your list using `[ ]` notation. [Difficulty: **Beginner**]

**solution:**
  ```{r}
mylist[4] 
mylist[[4]] 
```

39. Make a factor using factor(), with 5 elements.
Ex: `fa=factor(c("a","a","b"))`. [Difficulty: **Beginner**]

**solution:**
  ```{r}
fa=factor(c("a","a","b","c","c"))

```

40. Convert a character vector to a factor using `as.factor()`.
First, make a character vector using `c()` then use `as.factor()`.
[Difficulty: **Intermediate**]

**solution:**
  ```{r}
my.vec=c("a","a","b","c","c")
fa=as.factor(my.vec)
fa

```

41. Convert the factor you made above to a character using `as.character()`. [Difficulty: **Beginner**]

**solution:**
  ```{r}
fa
as.character(fa)

```



### Reading in and writing data out in R

1. Read CpG island (CpGi) data from the compGenomRData package `CpGi.table.hg18.txt`. This is a tab-separated file. Store it in a variable called `cpgi`. Use 
```
cpgFilePath=system.file("extdata",
                        "CpGi.table.hg18.txt",
                        package="compGenomRData")
```
to get the file path within the installed `compGenomRData` package. [Difficulty: **Beginner**]

**solution:**
  ```{r}
cpgFilePath
cpgi=read.table(file=cpgFilePath,header=TRUE,sep="\t")

```

2. Use `head()` on CpGi to see the first few rows. [Difficulty: **Beginner**]

**solution:**
  ```{r}
head(cpgi)

```

3. Why doesn't the following work? See `sep` argument at `help(read.table)`. [Difficulty: **Beginner**]

```{r readCpGex, eval=FALSE}
cpgtFilePath=system.file("extdata",
                "CpGi.table.hg18.txt",
                package="compGenomRData")
cpgtFilePath
cpgiSepComma=read.table(cpgtFilePath,header=TRUE,sep=",")
head(cpgiSepComma)
```

**solution:**
```{r}
cpgiSepComma=read.table("../data/CpGi.table.hg18.txt",header=TRUE,sep=",")
head(cpgiSepComma) 
```

4. What happens when you set `stringsAsFactors=FALSE` in `read.table()`? [Difficulty: **Beginner**]
```{r readCpGex12, eval=FALSE}
cpgiHF=read.table("intro2R_data/data/CpGi.table.hg18.txt",
                     header=FALSE,sep="\t",
                     stringsAsFactors=FALSE)
```                

**solution:**
The character column is now read as character instead of factor. 
```{r}
head(cpgiHF)
head(cpgi)
class(cpgiHF$V2)
class(cpgiHF$V2)
 
```
               

5. Read only the first 10 rows of the CpGi table. [Difficulty: **Beginner/Intermediate**]  

**solution:**
```{r}
cpgi10row=read.table("../data/CpGi.table.hg18.txt",header=TRUE,sep="\t",nrow=10)
cpgi10row 
```

6. Use `cpgFilePath=system.file("extdata","CpGi.table.hg18.txt",`
`package="compGenomRData")` to get the file path, then use
`read.table()` with argument `header=FALSE`. Use `head()` to see the results. [Difficulty: **Beginner**]  


**solution:**
```{r}
df=read.table(cpgFilePath,header=FALSE,sep="\t")
head(df) 
```

7. Write CpG islands to a text file called "my.cpgi.file.txt". Write the file
to your home folder; you can use `file="~/my.cpgi.file.txt"` in linux. `~/` denotes
home folder.[Difficulty: **Beginner**]  

**solution:**
```{r}
write.table(cpgi,file="~/my.cpgi.file.txt")
 
```


8. Same as above but this time make sure to use the `quote=FALSE`,`sep="\t"` and `row.names=FALSE` arguments. Save the file to "my.cpgi.file2.txt" and compare it with "my.cpgi.file.txt". [Difficulty: **Beginner**]  

**solution:**
```{r}
write.table(cpgi,file="~/my.cpgi.file2.txt",quote=FALSE,sep="\t",row.names=FALSE)
 
```


9. Write out the first 10 rows of the `cpgi` data frame.
**HINT:** Use subsetting for data frames we learned before. [Difficulty: **Beginner**]  

**solution:**
```{r}
write.table(cpgi[1:10,],file="~/my.cpgi.fileNrow10.txt",quote=FALSE,sep="\t")
 
```



10. Write the first 3 columns of the `cpgi` data frame. [Difficulty: **Beginner**] 

**solution:**
```{r}
dfnew=cpgi[,1:3]
write.table(dfnew,file="~/my.cpgi.fileCol3.txt",quote=FALSE,sep="\t") 
```

11. Write CpG islands only on chr1. **HINT:** Use subsetting with `[]`, feed a logical vector using `==` operator.[Difficulty: **Beginner/Intermediate**] 

**solution:**
```{r}
write.table(cpgi[cpgi$chrom == "chr1",],file="~/my.cpgi.fileChr1.txt",
            quote=FALSE,sep="\t",row.names=FALSE)
head(cpgi[cpgi$chrom == "chr1",])

 
```


12. Read two other data sets "rn4.refseq.bed" and "rn4.refseq2name.txt" with `header=FALSE`, and assign them to df1 and df2 respectively.
They are again included in the compGenomRData package, and you
can use the `system.file()` function to get the file paths. [Difficulty: **Beginner**] 

**solution:**
```{r}
df1=read.table("../data/rn4.refseq.bed",sep="\t",header=FALSE)
df2=read.table("../data/rn4.refseq2name.txt",sep="\t",header=FALSE)
```


13. Use `head()` to see what is inside the data frames above. [Difficulty: **Beginner**] 

**solution:**
```{r}
head(df1)
head(df2) 
```

14. Merge data sets using `merge()` and assign the results to a variable named 'new.df', and use `head()` to see the results. [Difficulty: **Intermediate**] 

**solution:**
```{r}
new.df=merge(df1,df2,by.x="V4",by.y="V1")
head(new.df)
 
```



### Plotting in R


Please run the following code snippet for the rest of the exercises. 
```{r plotExSeed}
set.seed(1001)
x1=1:100+rnorm(100,mean=0,sd=15)
y1=1:100
```

1. Make a scatter plot using the `x1` and `y1` vectors generated above. [Difficulty: **Beginner**] 

**solution:**
```{r}
plot(x1,y1)
 
```


2. Use the `main` argument to give a title to `plot()` as in `plot(x,y,main="title")`. [Difficulty: **Beginner**] 

**solution:**
```{r}
plot(x1,y1,main="scatter plot")
 
```


3. Use the `xlab` argument to set a label for the x-axis. Use `ylab` argument to set a label for the y-axis. [Difficulty: **Beginner**] 

**solution:**
```{r}
plot(x1,y1,main="scatter plot",xlab="x label")
 
```

4. Once you have the plot, run the following expression in R console. `mtext(side=3,text="hi there")` does. **HINT:** `mtext` stands for margin text. [Difficulty: **Beginner**] 

**solution:**
```{r}
plot(x1,y1,main="scatter plot",xlab="x label",ylab="y label")
mtext(side=3,text="hi there") 
 
```

5. See what `mtext(side=2,text="hi there")` does. Check your plot after execution. [Difficulty: **Beginner**] 

**solution:**
```{r}
mtext(side=2,text="hi there") 
mtext(side=4,text="hi there") 
 
```

6. Use *mtext()*  and *paste()* to put a margin text on the plot. You can use `paste()` as 'text' argument in `mtext()`. **HINT:** `mtext(side=3,text=paste(...))`. See how `paste()` is used for below. [Difficulty: **Beginner/Intermediate**] 

```{r pasteExample}
paste("Text","here")
myText=paste("Text","here")
myText
```

**solution:**
```{r}
mtext(side=3,text=paste("here","here"))

 
```


7. `cor()` calculates the correlation between two vectors.
Pearson correlation is a measure of the linear correlation (dependence) 
between two variables X and Y. Try using the `cor()` function on the `x1` and `y1` variables. [Difficulty: **Intermediate**] 

**solution:**
```{r}
corxy=cor(x1,y1) # calculates pearson correlation
 
```

8. Try to use `mtext()`,`cor()` and `paste()` to display the correlation coefficient on your scatter plot.  [Difficulty: **Intermediate**] 

**solution:**
```{r}
plot(x1,y1,main="scatter")
corxy=cor(x1,y1) 
#mtext(side=3,text=paste("Pearson Corr.",corxy))
mtext(side=3,text=paste("Pearson Corr.",round(corxy,3) ) )

plot(x1,y1)
mtext(side=3,text=paste("Pearson Corr.",round( cor(x1,y1)  ,3)  ) ) 
```

9. Change the colors of your plot using the `col` argument.
Ex: `plot(x,y,col="red")`. [Difficulty: **Beginner**] 

**solution:**
```{r}
plot(x1,y1,col="red")
 
```

10. Use `pch=19` as an argument in your `plot()` command. [Difficulty: **Beginner**] 

**solution:**
```{r}
plot(x1,y1,col="red",pch=19)
 
```


11. Use `pch=18` as an argument to your `plot()` command. [Difficulty: **Beginner**] 

**solution:**
```{r}
plot(x1,y1,col="red",pch=18)
?points 
```

12. Make a histogram of `x1` with  the `hist()` function. A histogram is a graphical representation of the data distribution. [Difficulty: **Beginner**] 

**solution:**
```{r}
hist(x1)

```


13. You can change colors with 'col', add labels with 'xlab', 'ylab', and add a 'title' with 'main' arguments. Try all these in a histogram.
[Difficulty: **Beginner**]  

**solution:**
```{r}
hist(x1, col = "red", xlab = "Distribution of X1", ylab = "Frequency Distribution", main = "Histogram of X1")
 
```

14. Make a boxplot of y1 with `boxplot()`.[Difficulty: **Beginner**] 

**solution:**
```{r}
boxplot(y1,main="title")
 
```

15. Make boxplots of `x1` and `y1` vectors in the same plot.[Difficulty: **Beginner**] 

**solution:**
```{r}
boxplot(x1,y1,ylab="values",main="title")
 
```

16. In boxplot, use the `horizontal = TRUE`  argument. [Difficulty: **Beginner**] 

**solution:**
```{r}
boxplot(x1,y1,ylab="values",main="title",horizontal=TRUE)
 
```

17. Make multiple plots with `par(mfrow=c(2,1))`
    - run `par(mfrow=c(2,1))` 
    - make a boxplot 
    - make a histogram
[Difficulty: **Beginner/Intermediate**] 

**solution:**
```{r}
par( mfrow=c(1,2) )
hist(x1)
boxplot(y1) 
```

18. Do the same as above but this time with `par(mfrow=c(1,2))`. [Difficulty: **Beginner/Intermediate**] 

**solution:**
```{r}
par(mfrow=c(2,2))
hist(x1)
boxplot(y1) 
```

19. Save your plot using the "Export" button in Rstudio. [Difficulty: **Beginner**] 

**solution:**
find and press Export button
 

20. You can make a scatter plot showing the density
of points rather than points themselves. If you use points it looks like this:

```{r colorScatterEx,out.width='50%'}

x2=1:1000+rnorm(1000,mean=0,sd=200)
y2=1:1000
plot(x2,y2,pch=19,col="blue")
```

If you use the `smoothScatter()` function, you get the densities. 
```{r smoothScatterEx,out.width='50%'}
smoothScatter(x2,y2,
              colramp=colorRampPalette(c("white","blue",
                                         "green","yellow","red"))) 
```
 
Now, plot with the `colramp=heat.colors` argument and then use a custom color scale using the following argument.
```
colramp = colorRampPalette(c("white","blue", "green","yellow","red")))
```
[Difficulty: **Beginner/Intermediate**] 

**solution:**
```{r}
smoothScatter(x2,y2,colramp = heat.colors )
smoothScatter(x2,y2,
               colramp = colorRampPalette(c("white","blue", "green","yellow","red")))
 
```

### Functions and control structures (for, if/else, etc.)
Read CpG island data as shown below for the rest of the exercises.

```{r CpGexReadchp2,eval=TRUE}
cpgtFilePath=system.file("extdata",
                "CpGi.table.hg18.txt",
                package="compGenomRData")
cpgi=read.table(cpgtFilePath,header=TRUE,sep="\t")
head(cpgi)
```

1. Check values in the perGc column using a histogram.
The 'perGc' column in the data stands for GC percent => percentage of C+G nucleotides. [Difficulty: **Beginner**] 

**solution:**
```{r}
hist(cpgi$perGc) # most values are between 60 and 70
 
```

2. Make a boxplot for the 'perGc' column. [Difficulty: **Beginner**] 

**solution:**
```{r}
boxplot(cpgi$perGc) 
 
```


3. Use if/else structure to decide if the given GC percent is high, low or medium.
If it is low, high, or medium: low < 60, high>75, medium is between 60 and 75;
use greater or less than operators, `<`  or ` >`. Fill in the values in the code below, where it is written 'YOU_FILL_IN'. [Difficulty: **Intermediate**]
```{r functionEvExchp2}

GCper=65

  # check if GC value is lower than 60, 
  # assign "low" to result
  if('YOU_FILL_IN'){
    result="low"
    cat("low")
  }
  else if('YOU_FILL_IN'){  # check if GC value is higher than 75,      
                           #assign "high" to result
    result="high"
    cat("high")
  }else{ # if those two conditions fail then it must be "medium"
    result="medium"
  }

result

```

**solution:**
```{r}
GCper=65
 #result="low"# set initial value
  
  if(GCper < 60){ # check if GC value is lower than 60, assign "low" to result
    result="low"
    cat("low")
  }
  else if(GCper > 75){  # check if GC value is higher than 75, assign "high" to result
    result="high"
    cat("high")
  }else{ # if those two conditions fail then it must be "medium"
    result="medium"
  }

result

 
```

4. Write a function that takes a value of GC percent and decides
if it is low, high, or medium: low < 60, high>75, medium is between 60 and 75.
Fill in the values in the code below, where it is written 'YOU_FILL_IN'. [Difficulty: **Intermediate/Advanced**]


```
GCclass<-function(my.gc){
  
  YOU_FILL_IN
  
  return(result)
}
GCclass(10) # should return "low"
GCclass(90) # should return "high"
GCclass(65) # should return "medium"
```

**solution:**
```{r}
GCclass<-function(my.gc){
  
  result="low"# set initial value
  
  if(my.gc < 60){ # check if GC value is lower than 60, assign "low" to result
    result="low"
  }
  else if(my.gc > 75){  # check if GC value is higher than 75, assign "high" to result
    result="high"
  }else{ # if those two conditions fail then it must be "medium"
    result="medium"
  }
  return(result)
}
GCclass(10) # should return "low"
GCclass(90) # should return "high"
GCclass(65) # should return "medium" 
```


5. Use a for loop to get GC percentage classes for `gcValues` below. Use the function
you wrote above.[Difficulty: **Intermediate/Advanced**]

```
gcValues=c(10,50,70,65,90)
for( i in YOU_FILL_IN){
  YOU_FILL_IN
}
```

**solution:**
```{r}
gcValues=c(10,50,70,65,90)
for( i in gcValues){
 
  print(GCclass(i) )
} 
```


6. Use `lapply` to get GC percentage classes for `gcValues`. [Difficulty: **Intermediate/Advanced**]

```{r lapplyExExerciseChp2,eval=FALSE}
vec=c(1,2,4,5)
power2=function(x){ return(x^2)  }
    lapply(vec,power2)
```


**solution:**
```{r}
s=lapply(gcValues,GCclass)
 
```



7. Use sapply to get values to get GC percentage classes for `gcValues`. [Difficulty: **Intermediate**]

**solution:**
```{r}
s=sapply(gcValues,GCclass)
 
```

8. Is there a way to decide on the GC percentage class of a given vector of `GCpercentages`
without using if/else structure and loops ? if so, how can you do it?
**HINT:** Subsetting using < and > operators.
[Difficulty: **Intermediate**]

**solution:**
```{r}
result=rep("low",length(gcValues) )
result[gcValues > 75]="high"
result[gcValues < 75 & gcValues > 60 ] = "medium" 
```
