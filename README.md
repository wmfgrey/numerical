
A program for automatically generating school reports

This program writes school reports.  It requires a list of comments in a comments bank.  Also a set of grades (1-5) for various student aptutides are needed, where 1 is the highest grade and 5 the lowest grade for a given criteria.  If there is more than one comment per grade and criteria the code randomly chooses a suitable comment for each pupils. The maximum number of criteria that you can have is 13. Example datasets are given in the folder dataFiles. The program needs to be run twice. Firstly the user needs to give pupils mark for the chosen criteria given in the comments bank data file.  This can include effort, behaviour, attainment and so on depending on how you choose to set it up. Secondly you can generate reports for each pupil based on marks given.  For instance, for the first run you might use: 

	reportGenerator -grades classList.txt bank.xml classListOut.txt 

This requires you to have a text file (e.g. classList.txt) with a list of students needs to be in the format:

<firstname> <lastname> <Gender M or F>
For instance,

	Bart Simpson M 
	Lisa Simpson F 
	Nelson Muntz M
	Ralf Wiggum M 

YOu will also need the comments bank file.  There is an example in the dataFiles called bank.xml.  You will need to edit this to include your own comments.

An output class list file will be generated with grades and should look similar to (e.g. classListOut.txt):

	Bart Simpson M 1 1 1 1 1 1 1 1 1 1 1 1
	Lisa Simpson F 1 1 1 1 1 1 1 1 1 1 1 1
	Nelson Muntz M 2 2 2 2 2 2 2 2 2 2 2 2
	Ralf Wiggum M 5 5 5 5 5 5 5 5 5 5 5 5

As the program is run you will need to enter grades for each pupil according to each criteria.  This can be quite time consuming so you could edit the file directly, but you will need to know which column corresponds to which criteria. 

On the second run of the program you will run:
	
	reportGenerator -reports classListOut.txt bank.xml reports.txt

The final output file is reports.txt and should include reports for all pupils in the original class list. The writeReport wrapper script will help you to get started.





