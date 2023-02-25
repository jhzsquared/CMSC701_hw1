// java -cp ":/home/jzhu/umd/CMSC701/hw1_suffixarray/lib/protobuf-java-3.21.12.jar" buildsa "/home/jzhu/umd/CMSC701/hw1_suffixarray/data/ecoli.fa" "/home/jzhu/umd/CMSC701/hw1_suffixarray/data/ecoli_sa.bin"
// javac -cp ":/home/jzhu/umd/CMSC701/hw1_suffixarray/lib/protobuf-java-3.21.12.jar" buildsa.java

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

public class buildsa {
    
    static String ReadFastaFile(String reference){
         //flie path needs to be full path   
         // assumes onee sequence per file  
        boolean firstLine = true;
        String genome;
        StringBuilder stringBuilder = new StringBuilder();

        try (BufferedReader reader = new BufferedReader(new FileReader(reference))){
            String line = reader.readLine();
            while (line != null) {
                if (line.charAt(0) == '>' && firstLine == true) {
                    stringBuilder = new StringBuilder();
                    firstLine = false;
                } else if (line.charAt(0) != '>'){
                    stringBuilder.append(line.toUpperCase());
                }
                line = reader.readLine();
            }
            reader.close();

        } catch (FileNotFoundException e){
            System.err.println(e);
            System.exit(1);
         } catch (IOException e)  {
			System.out.println(e);
			System.exit(1);
		}
       
        // Replace any Ns with random letter for sake of project
        int nIndex = stringBuilder.indexOf("N");
        while (nIndex > 0) {
            Random rand = new Random();
            char[] dna = {'A','C','G','T'};
            stringBuilder.setCharAt(nIndex, dna[(rand.nextInt(4))]);
            nIndex =  stringBuilder.indexOf("N", nIndex+1);
        }
        genome = stringBuilder.toString();
        return genome;
        
    }
    
    static public Map<String, SuffixArrayProto.indexInterval> CreatePrefixTable(int[] suffixArray, String genome, int k){
        // Create prefix table: key: prefix, value: [start-inclusive, end-exclusive]) 
        Map<String,SuffixArrayProto.indexInterval> prefixTable = new HashMap<String, SuffixArrayProto.indexInterval>(); 
        //initialize 
        String prefix = ""; 
        String new_prefix;
        int startIndex = 1;
        int i=1;
        boolean first = true;
        while (i < suffixArray.length){
            if (suffixArray[i]<=(genome.length()-k)){ //not too short substring
                new_prefix = genome.substring(suffixArray[i], suffixArray[i]+k);
                if (prefix.equals(new_prefix)){
                    i+=1; //Same prefix
                } else{ //different prefix, store previous prefix interval into hashmap
                    if (first){//don't save
                    } else{ //save results 
                        if (prefix.length()==k){ //so we can ignore the short prefixes
                            prefixTable.put(prefix, SuffixArrayProto.indexInterval.newBuilder().addAllIndices(new ArrayList<>(Arrays.asList(startIndex, i))).build());
                        }
                    }
                    //update values
                    startIndex = i;
                    prefix = new_prefix;
                    i+=1;
                    first = false;
                }          
            } else{ //save results and continue to next suffix (use short one as a placeholder)
                if (prefix.length()==k){
                    prefixTable.put(prefix, SuffixArrayProto.indexInterval.newBuilder().addAllIndices(new ArrayList<>(Arrays.asList(startIndex, i))).build());   
                } else {
                    //Don't save just move on
                } 
                first = false;
                prefix = genome.substring(suffixArray[i]);
                startIndex = i;
                i+=1; //skip to next suffix
            }         
        }
        // save last prefix;     
        prefixTable.put(prefix, SuffixArrayProto.indexInterval.newBuilder().addAllIndices(new ArrayList<>(Arrays.asList(startIndex, i))).build());    
        return prefixTable;
    }

    public static void main(String[] args) {
        int k=0;
        String reference;
        String output;
        boolean preftab = false;
        if (args.length == 2) {
            reference = args[0];
            output = args[1];
        } else if (args.length == 4) {//preftab request (args[0] should have been --preftab)
            preftab = true;
            k = Integer.parseInt(args[1]);
            reference = args[2];
            output = args[3];
        } else {
            System.err.println("Invalid arguments count: " + args.length);
            return;
        }
        // read in fasta file
        String genome = ReadFastaFile(reference);
        //add "$" if it doesn't start with it
        if (genome.charAt(genome.length()-1) != "$".charAt(0)) {
            genome = genome+"$";
        }

        //set up for suffix array construction
        char[] chars = genome.toCharArray();
        int[] s = new int[chars.length];
        int minimum = 256;
        for (int j = 0; j< chars.length; j++) {
            s[j] = chars[j] + 1;
            if (s[j] < minimum) 
                minimum = s[j];
        }
        int maximum = 0;
        for (int j = 0; j< chars.length; j++) {
            s[j] = chars[j] + 2 - minimum;		
            if (s[j] > maximum)
                maximum = s[j];
        }
        
        // time suffix array construction
        long startTime = System.nanoTime();
        final int[] suffixArray = SuffixArray.constructSuffixArray(s, maximum);	
        long endTime = System.nanoTime();
        long duration = (endTime-startTime);
        System.out.println("Constructing the suffix array took "+ duration);
        
        //Encode results to binary via protobuf
        final SuffixArrayProto.SuffixMessage.Builder messageBuilder = SuffixArrayProto.SuffixMessage.newBuilder();
        messageBuilder.setGenome(genome);
        //convert SA to iterable so it can be saved
        Iterable<Integer> iterSuffix = () -> IntStream.of(suffixArray).boxed().iterator();
        messageBuilder.addAllSuffixArray(iterSuffix);

        if (preftab){
            //create prefix table
            long prefixStartTime = System.nanoTime();
            Map<String,SuffixArrayProto.indexInterval> prefixTable = CreatePrefixTable(suffixArray, genome, k);	
            long prefixEndTime = System.nanoTime();
            long prefixDuration = (prefixEndTime-prefixStartTime);
            System.out.println("Constructing the prefix table took "+ prefixDuration);
            // //save to protobuf
            messageBuilder.putAllPreftab(prefixTable);
        }
        //save protobuf binary file
        try(FileOutputStream outputstream = new FileOutputStream(output)){
            messageBuilder.build().writeTo(outputstream);
            outputstream.close();    
        } catch (IOException e) {
            System.out.println(e);
        }
        // // test if it's right by reading and comparing
        // try (FileInputStream input = new FileInputStream(output)) {
        //     SuffixArrayProto.SuffixMessage.Builder messageRead = SuffixArrayProto.SuffixMessage.newBuilder();
        //     SuffixArrayProto.SuffixMessage mess = messageRead.mergeFrom(input).build();
        //     System.out.println(mess.getSuffixArray(1));
        // } catch (IOException e) {
        //     System.out.println(e);
        // }

    }
}   