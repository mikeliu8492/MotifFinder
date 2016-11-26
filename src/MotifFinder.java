import java.util.*;
import java.io.*;


public class MotifFinder
{

	/**
	 * Global variables
	 * 
	 */
	
	static int MOTIF_LENGTH = 0;
	static String SEQUENCE_FILE = "sequences.fa";
	static String LENGTH_FILE = "motiflength.txt";
	
	//list of best motif start locations, 0-indexed
	static ArrayList<Integer> bestMotifLocations = new ArrayList<Integer>();
	
	
	//list of DNA sequences of length 500
	static ArrayList<String> sequenceList = new ArrayList<String>();
	
	/**
	 * This reads in the sequences in the "file" and adds them to the "sequenceList" ArrayList.
	 * 
	 * @param file		file of sequences to read in
	 * @throws IOException
	 */
	static void formSequenceList(String file) throws IOException
	{
		int idx = 0;
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		 
		String line = null;
		while ((line = br.readLine()) != null)
		{
			if(line.length() == 0)
				continue;
			else if (line.charAt(0) == '>')
				continue;
			else
				sequenceList.add(line);
		}
	
	}
	
	
	
	
	
	
	/**
	 * function to write out locations of the files Array +1
	 * @throws IOException
	 */
	public static void writeOutSites() throws IOException
	{
		int listSize = bestMotifLocations.size();
		
		FileWriter write = new FileWriter("predictedsites.txt", false);
		PrintWriter print_line = new PrintWriter(write);
		
		for(int i = 0; i < listSize; i++)
		{
			print_line.println(">Sequence"+i);
			int trueIdx = bestMotifLocations.get(i) + 1;
			print_line.println(trueIdx);
			print_line.println();
		}
		
		print_line.close();
	}
	
	
	/**
	 * Write out the tab-delimited motif file similar to part 1.
	 * 
	 * @throws IOException
	 */
	public static void writeOutMotif() throws IOException
	{
		int listSize = bestMotifLocations.size();
		
		FileWriter write = new FileWriter("predictedmotifs.txt", false);
		PrintWriter print_line = new PrintWriter(write);
		
		int rows = MOTIF_LENGTH;
		int col = 4; 
		int matrix[][] = new int[rows][4];
		
		for (int m = 0; m < listSize; m++)
		{
			int start = bestMotifLocations.get(m);
			String currentMotif = sequenceList.get(m).substring(start, start+MOTIF_LENGTH);
			
			for(int i = 0; i < rows; i++)
			{
				if(currentMotif.charAt(i) == 'A')
					matrix[i][0] += 1; 
				else if(currentMotif.charAt(i) == 'C')
					matrix[i][1] += 1; 
				else if(currentMotif.charAt(i) == 'G')
					matrix[i][2] += 1; 
				else
					matrix[i][3] += 1; 
			}
		}
		
		print_line.println(">MOTIF1 " + MOTIF_LENGTH);
		
		for(int i = 0; i < rows; i++)
		{
			print_line.printf("%d\t%d\t%d\t%d\n", matrix[i][0], matrix[i][1], matrix[i][2], matrix[i][3]);
		}
		
		
		print_line.println("<");
		
		
		print_line.close();
	}
	
	

	/**
	 * Takes the motif-length file and parses out the motif-length
	 * 
	 * @throws IOException
	 */
	public static void importLengthFromFile(String file) throws IOException
	{
		int idx = 0;
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		 
		String line = null;
		line = br.readLine();
		
		MOTIF_LENGTH = Integer.parseInt(line);
	}

	/**
	 * This takes the minimum-hamming distance between all motifs in a given ArrayList.
	 * It looks at each position of motifs and calculates the most common base (e.g. A, C, G, T).
	 * Any base that deviates from that most common base ina given positionwill be added to the 
	 * overall hamming-distances core.  The larger the score, the less consensus there is between sequences.
	 * 
	 * 
	 * 
	 * @paramArray  this represents the array of motif start indices
	 *
	 * @return 	int representing the minimum Hamming distance score
	 */
	public static int getHammingScore(ArrayList<Integer> paramArray)
	{
		
		ArrayList<Integer> evaluateArray = new ArrayList<Integer>();
		for(Integer obj : paramArray)
		{
		    evaluateArray.add(obj);
		}
		
		int [][] table = new int [MOTIF_LENGTH][4];
		int totalScore = 0;
		

		for (int i = 0; i < evaluateArray.size(); i++)
		{
			int startLocation = evaluateArray.get(i);

			String myString = sequenceList.get(i).substring(startLocation, startLocation+MOTIF_LENGTH);

			for(int j = 0; j <MOTIF_LENGTH;j++)
			{
				if(myString.charAt(j) == 'A')
					table[j][0]++;
				else if(myString.charAt(j) == 'C')
					table[j][1]++;
				else if(myString.charAt(j) == 'G')
					table[j][2]++;
				else
					table[j][3]++;

			}
			
		}

		

		for (int i = 0; i < MOTIF_LENGTH; i++)
		{
			int max = -1;
			int idx = 0;
			int localTotal = 0;

			for(int m = 0; m < 4; m++)
			{
				if(table[i][m] > max)
				{
					max = table[i][m];
					idx = m;
				}
					
			}

			for(int m = 0; m < 4; m++)
			{
				if(m != idx)
					localTotal += table[i][m];
			}
			
			totalScore += localTotal;
		}

		return totalScore;
	}
	
	public static void main(String[] args) throws IOException
	{
		formSequenceList(args[0]);
		importLengthFromFile(args[1]);
		
		
		/*
		 Take the 0-th index of each sequence and enter the motif.
		*/
		for(int i = 0; i < sequenceList.size(); i++)
			bestMotifLocations.add(0);

		//set baseline minimum score
		int minScore = getHammingScore(bestMotifLocations);
		int maxIndex = sequenceList.get(0).length()-MOTIF_LENGTH+1;
		

		for(int firstRowIdx = 0; firstRowIdx < maxIndex; firstRowIdx++)
		{
			//create a local working copy of an array-list of motifs
			ArrayList<Integer> local = new ArrayList<Integer>();
			local.add(firstRowIdx);
			
			/*
			for each other sequence in list, select the substring that will minimize overall hamming distance score.
			keep doing this to assemble an ArrayList of locations that will resemble the most likely locations of 
			the motifs and have the best consensus
			*/

			for(int otherRows = 1; otherRows < sequenceList.size(); otherRows++)
			{
				int stepMin = Integer.MAX_VALUE;
				int targetIdx = 0;
				

				for(int otherRowIdx = 0; otherRowIdx < maxIndex; otherRowIdx++)
				{
					
					local.add(otherRowIdx);
					
					int scoreLocal = getHammingScore(local);
					if(scoreLocal < stepMin)
					{
						stepMin = scoreLocal;
						targetIdx = otherRowIdx;
					}
					
					int lastIdx = local.size()-1;
					local.remove(lastIdx);
				}
				
				local.add(targetIdx);
			}
			/*
			if score of your new list of best motif location is better and old score
			replace the minScore variable you set with the new score
			 and replace the bestMotifLocations list with your new list
			*/
			
			int localScore = getHammingScore(local);
			if (localScore < minScore)
			{
				minScore = localScore;
				bestMotifLocations = local;
			}
		}
		

		
		for(int i = 0; i < bestMotifLocations.size(); i++)
		{
			int location = bestMotifLocations.get(i)+1;
			System.out.println("Location:  " + location);
		}
		
		
		writeOutSites();
		writeOutMotif();
		
	}
}
