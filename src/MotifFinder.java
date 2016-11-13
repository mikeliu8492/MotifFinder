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
	
	//list of motif start locations 0-indexed
	static ArrayList<Integer> motifLocations = new ArrayList<Integer>();
	
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
	 * 
	 * Score the PWN for the first two sequences.
	 * 
	 * @param first		substring from the first DNA sequence
	 * @param second	substring from the second DNA sequence
	 * 
	 * We assume a background rate of all 0.25, since all bases as just as likely to appear. 
	 * 
	 * @return
	 */
	static double scorePWN(String first, String second)
	{
		double pwn = 0;
		
		int localChart [][] = new int[4][MOTIF_LENGTH];
		
		for (int i = 0; i < first.length(); i++)
		{
			if (first.charAt(i) == 'A')
				++localChart[0][i];
			else if (first.charAt(i) == 'C')
				++localChart[1][i];
			else if (first.charAt(i) == 'G')
				++localChart[2][i];
			else
				++localChart[3][i];
			
			if (second.charAt(i) == 'A')
				++localChart[0][i];
			else if (second.charAt(i) == 'C')
				++localChart[1][i];
			else if (second.charAt(i) == 'G')
				++localChart[2][i];
			else
				++localChart[3][i];
			
		}
		
		for (int i = 0; i < MOTIF_LENGTH; i++)
		{
			double local_sum = 0;
			for (int m = 0; m < 4; m++)
			{
				double p = (double)(localChart[m][i])/4;
				if (p != 0)
				{
					double result = p;
					double division = Math.log(p/0.25)/Math.log(2);
					result *= division;
					local_sum += result;
				}
			}
			
			pwn += local_sum;
			
		}
		
		
		return pwn;
	}
	
	
	
	/**
	 * Scoring function for all the others
	 * @param target
	 * @return
	 */
	static double scoreAllOthers()
	{
		int width = MOTIF_LENGTH;
		int rows = 4; 
		int matrix[][] = new int[rows][width];
		
		for (int m = 0; m < motifLocations.size(); m++)
		{
			int start = motifLocations.get(m);
			String currentMotif = sequenceList.get(m).substring(start, start+MOTIF_LENGTH);
			
			for(int i = 0; i < width; i++)
			{
				if(currentMotif.charAt(i) == 'A')
					matrix[0][i]++; 
				else if(currentMotif.charAt(i) == 'C')
					matrix[1][i]++; 
				else if(currentMotif.charAt(i) == 'G')
					matrix[2][i]++; 
				else
					matrix[3][i]++; 
			}
		}

		
		double final_score = 0;
		for (int i = 0; i < width; i++)
		{
			final_score += Math.max(Math.max(matrix[0][i], matrix[1][i]) , Math.max(matrix[2][i], matrix[3][i])); 
		}
		
		
		return final_score; 
		
	}
	
	
	/**
	 * function to write out locations of the files Array +1
	 * @throws IOException
	 */
	public static void writeOutSites() throws IOException
	{
		int listSize = motifLocations.size();
		
		FileWriter write = new FileWriter("predictedsites.txt", false);
		PrintWriter print_line = new PrintWriter(write);
		
		for(int i = 0; i < listSize; i++)
		{
			print_line.println(">Sequence"+i);
			int trueIdx = motifLocations.get(i) + 1;
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
		int listSize = motifLocations.size();
		
		FileWriter write = new FileWriter("predictedmotifs.txt", false);
		PrintWriter print_line = new PrintWriter(write);
		
		int width = MOTIF_LENGTH;
		int length = 4; 
		int matrix[][] = new int[length][width];
		
		for (int m = 0; m < listSize; m++)
		{
			int start = motifLocations.get(m);
			String currentMotif = sequenceList.get(m).substring(start, start+MOTIF_LENGTH);
			
			for(int i = 0; i < width; i++)
			{
				if(currentMotif.charAt(i) == 'A')
					matrix[0][i] += 1; 
				else if(currentMotif.charAt(i) == 'C')
					matrix[1][i] += 1; 
				else if(currentMotif.charAt(i) == 'G')
					matrix[2][i] += 1; 
				else
					matrix[3][i] += 1; 
			}
		}
		
		print_line.println(">MOTIF1 " + MOTIF_LENGTH);
		
		for(int i = 0; i < width; i++)
		{
			print_line.printf("%d\t%d\t%d\t%d\n", matrix[0][i], matrix[1][i], matrix[2][i], matrix[3][i]);
		}
		
		
		print_line.println("<");
		
		
		print_line.close();
	}
	
	
	public static void importLengthFromFile(String file) throws IOException
	{
		int idx = 0;
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		 
		String line = null;
		line = br.readLine();
		
		MOTIF_LENGTH = Integer.parseInt(line);
	}
	
	public static void main(String[] args) throws IOException
	{
		formSequenceList(SEQUENCE_FILE);
		importLengthFromFile(LENGTH_FILE);
		
		String first = sequenceList.get(0);
		String second = sequenceList.get(1);
		
		
		String sub1 = null;
		String sub2 = null;
		
		double best = -1;
		int firstIdx = 0;
		int secondIdx = 0;
		
		
		//find best match between 1st and 2nd substring and construct a matrix
		
		for (int i = 0; i < first.length() - MOTIF_LENGTH; i++)
		{
			sub1 = first.substring(i, i+MOTIF_LENGTH);
			for(int m = 0; m < second.length() - MOTIF_LENGTH; m++)
			{
				sub2 = second.substring(m, m+MOTIF_LENGTH);
				
				double potential = scorePWN(sub1, sub2);
				
				if(potential >= best)
				{
					best = potential;
					firstIdx = i;
					secondIdx = m;
				}
					
			}
		}		
		
		if(firstIdx > 0)
			firstIdx -= 1;
		
		if(secondIdx > 0)
			secondIdx -= 1;
		
		motifLocations.add(firstIdx);
		motifLocations.add(secondIdx);
		
		
			
		for(int i = 2; i < sequenceList.size(); i++)
		{
			
			String current = sequenceList.get(i);
			
			int local_idx = 0;
			double best_score = -1;
			
			for (int m = 0; m < current.length()-MOTIF_LENGTH+1; m++)
			{
				//add the current location of substring start to calculation matrix
				motifLocations.add(m);
				double currentScore = scoreAllOthers();
				if (currentScore >= best_score)
				{
					best_score = currentScore;
					local_idx = m;
				}
				
				//remove from calculation matrix
				motifLocations.remove(motifLocations.size()-1);
						
			}
			
			//add the ideal substring to calculation matrix
			motifLocations.add(local_idx);
		}
		
		
		/*
		for (int m = 0; m < motifLocations.size(); m++)
		{
			int start = motifLocations.get(m);
			String currentMotif = sequenceList.get(m).substring(start-1, start+MOTIF_LENGTH);
			System.out.println(currentMotif);
		}
		*/
		
		writeOutSites();
		writeOutMotif();
		
	}
}
