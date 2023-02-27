//java -cp ":/home/jzhu/umd/CMSC701/hw1_suffixarray/lib/protobuf-java-3.21.12.jar" querysa /home/jzhu/umd/CMSC701/hw1_suffixarray/data/ecoli_sa.bin /home/jzhu/umd/CMSC701/hw1_suffixarray/data/fake_query.fa naive /home/jzhu/umd/CMSC701/hw1_suffixarray/data/query_output.bin
import java.io.*;
import java.util.*;
import java.util.Collections;
import java.util.stream.Collectors;

public class querysa {
    // read in binary string and suffix array (and maybe prefix table) 
    static SuffixArrayProto.SuffixMessage ReadSABin(String binPath){
        // Read in protobuf output
        SuffixArrayProto.SuffixMessage mess=null;
        try (FileInputStream input = new FileInputStream(binPath)) {
            SuffixArrayProto.SuffixMessage.Builder messageRead = SuffixArrayProto.SuffixMessage.newBuilder();
            mess = messageRead.mergeFrom(input).build();

        } catch (IOException e) {
            System.out.println(e);
            System.exit(1);
        }
        return mess;
    }

    static Map<String,String> ReadFastaQuery(String queryPath){
         //flie path needs to be full path   
         // assumes onee sequence per file 
        String query,header = "";
        boolean firstLine = true;
        StringBuilder stringBuilder = new StringBuilder();
        Map<String,String> queries = new HashMap<String, String>();

         try (BufferedReader reader = new BufferedReader(new FileReader(queryPath))){
             String line = reader.readLine();
             while (line != null) {
                 if (line.charAt(0) == '>' && firstLine) {
                    //new query coming up
                     stringBuilder = new StringBuilder();
                     String[] splitHeader = line.split("\\s+");
                     header = splitHeader[0].substring(1); // remove the >
                     firstLine = false;
                 } else if (line.charAt(0) != '>'){
                    //getting query lines
                     stringBuilder.append(line.toUpperCase());
                 } else if (line.charAt(0) == '>' && firstLine==false){
                    //new query request, so save previous reesults
                    query = stringBuilder.toString();
                    queries.put(header,query);
                    //initiate for next query
                    stringBuilder = new StringBuilder();
                    String[] splitHeader = line.split("\\s+");
                    header = splitHeader[0].substring(1); // remove the >
                 }
                 line = reader.readLine();
             }
             reader.close();
             //save final query result
            query = stringBuilder.toString();
            queries.put(header,query);
 
         } catch (FileNotFoundException e){
             System.err.println(e);
             System.exit(1);
          } catch (IOException e)  {
             System.out.println(e);
             System.exit(1);
         }
         return queries;
    }

    static String getSubstring(String genome, int start, int length){
        String substring;
        try{
            substring = genome.substring(start, start + length);
        } catch (java.lang.StringIndexOutOfBoundsException e){ //start is too short
            substring = genome.substring(start);
        }
        return substring;
    }
    

    static List<Integer> naiveQuery(String query, String genome, List<Integer> suffixArray, Map<String,SuffixArrayProto.indexInterval> prefixTable){
        // naive binary search (if prefix table is not null then it will use it)
        List<Integer> hits = new ArrayList<Integer>();
        int l= 0;
        int r = suffixArray.size()-1;
        int queryLength = query.length();
        int c = (int)Math.floor((l+r)/2);
        //start timing
        long startTime = System.nanoTime();

        if (prefixTable.size()!=0){
        // set initial left and right limits based on prefixtable findings
            //check if pattern is longer or shorter than prefixes
            int prefixLength = prefixTable.keySet().iterator().next().length();
            if (prefixLength > queryLength){
                //query is shorter than prefixes generated
                //find range of prefixes that fit //binary search --and then you have effectively found range
                //binary search for lower bound
                List<String> prefixList = new ArrayList<String>(prefixTable.keySet());
                Collections.sort(prefixList); //lexicographically sort
                l = 0;
                r = prefixList.size();
                while (l< r){
                    c = (int)Math.floor((l+r)/2);
                    //check middle
                    if (getSubstring(prefixList.get(c), 0, queryLength).compareTo(query) >= 0){
                        // query is  smaller or equal  to center substring
                        r = c; //left bisect
                    } else{
                        //query is bigger than center
                        l = c+1;// right bisect
                    }
                }
                int leftBound = l;
                //binary search for upper bound (might as well start at the left bound)
                r = prefixList.size();
                while (l<r){
                    c = (int)Math.floor((l+r)/2);
                    //check middle
                    if (getSubstring(prefixList.get(c), 0, queryLength).compareTo(query) <= 0){
                        // query is bigger or equal  to center substring
                        l = c+1; //right bisect
                    } else{
                        //query is bigger than center
                        r = c;// right bisect
                    }
                }
                // convert prefixtable hits to actual genome indices
                leftBound = prefixTable.get(prefixList.get(leftBound)).getIndices(0);
                r = prefixTable.get(prefixList.get(r)).getIndices(1);
                for (int i = leftBound; i<r; i++){
                    hits.add(suffixArray.get(i).intValue());
                }   
                //stop timing
                long endTime = System.nanoTime();
                long duration = (endTime-startTime);
                System.out.println("Query of "+ query.length() + " length took " + duration);
                return hits;
                
            } else { 
                //query is longer than prefix, proceed with standard search
                String queryPrefix = query.substring(0, prefixLength);
                if (prefixTable.get(queryPrefix)==null){
                    long endTime = System.nanoTime();
                    long duration = (endTime-startTime);
                    System.out.println("Query of "+ query.length() + " length took " + duration);
                    return hits; //no matches if prefix can not be found
                } else{
                    l = prefixTable.get(queryPrefix).getIndices(0);
                    r = prefixTable.get(queryPrefix).getIndices(1);
                }             
            }
        }
        //binary search for lower bound
        while (l< r){
            c = (int)Math.floor((l+r)/2);
            //check middle
            if (getSubstring(genome, suffixArray.get(c).intValue(),queryLength).compareTo(query) >= 0){
                // query is  smaller or equal  to center substring
                r = c; //left bisect
            } else{
                //query is bigger than center
                l = c+1;// right bisect
            }
        }
        int leftBound = l;
        //binary search for upper bound (might as well start at the left bound)
        r = suffixArray.size()-1;
        while (l<r){
            c = (int)Math.floor((l+r)/2);
            //check middle
            if (getSubstring(genome, suffixArray.get(c).intValue(),queryLength).compareTo(query) <= 0){
                // query is bigger or equal  to center substring
                l = c+1; //right bisect
            } else{
                //query is bigger than center
                r = c;// right bisect
            }
        }
        // convert suffix Array hits to actual genome indices
        for (int i = leftBound; i<r; i++){
            hits.add(suffixArray.get(i).intValue());
        }   
        //stop timing
        long endTime = System.nanoTime();
        long duration = (endTime-startTime);
        System.out.println("Query of "+ query.length() + " length took " + duration);
        return hits;
    }

    static int getLCP(String a, String b){
        //get length of longest common prefix between two strings
        int minSize = Math.min(a.length(), b.length());
        int i=0;
        while (i<minSize && a.charAt(i)==b.charAt(i)){
            i++;
        }
        String lcp = a.substring(0,i);
        return lcp.length();
    }

    static List<Integer> simpAccelQuery(String query, String genome, List<Integer> suffixArray, Map<String,SuffixArrayProto.indexInterval> prefixTable){
        //use binary search with LCPs to skip some character comparisons
        List<Integer> hits = new ArrayList<Integer>();
        int l= 0;
        int r = suffixArray.size()-1;
        int queryLength = query.length();
        int c = (int)Math.floor((l+r)/2);
        //start timing
        long startTime = System.nanoTime();

        if (prefixTable.size()!=0){
        // set initial left and right limits based on prefixtable findings
            //check if pattern is longer or shorter than prefixes
            int prefixLength = prefixTable.keySet().iterator().next().length();
            if (prefixLength > queryLength){
                //query is shorter than prefixes generated
                //find range of prefixes that fit //binary search --and then you have effectively found range
                //binary search for lower bound
                List<String> prefixList = new ArrayList<String>(prefixTable.keySet());
                Collections.sort(prefixList); //lexicographically sort
                l = 0;
                r = prefixList.size()-1;
                int lLcp = getLCP(getSubstring(prefixList.get(l), 0, queryLength), query);
                int rLcp = getLCP(getSubstring(prefixList.get(r), 0, queryLength), query);
                while (l< r){
                    c = (int)Math.floor((l+r)/2);
                    int minLcp = Math.min(lLcp, rLcp);
                    //check middle
                    if (getSubstring(prefixList.get(c), minLcp, queryLength-minLcp).compareTo(query.substring(minLcp)) >= 0){
                        // query is  smaller or equal  to center substring
                        r = c; //left bisect
                        rLcp = getLCP(getSubstring(prefixList.get(r), 0, queryLength), query);
                    } else{
                        //query is bigger than center
                        l = c+1;// right bisect
                        lLcp = getLCP(getSubstring(prefixList.get(l), 0, queryLength), query);
                    }
                }
                int leftBound = l;
                //binary search for upper bound (might as well start at the left bound)
                r = prefixList.size()-1;
                rLcp = getLCP(getSubstring(prefixList.get(r), 0, queryLength), query);
                while (l<r){
                    c = (int)Math.floor((l+r)/2); 
                    int minLcp = Math.min(lLcp, rLcp);               
                    //check middle
                    if (getSubstring(prefixList.get(c), minLcp, queryLength-minLcp).compareTo(query.substring(minLcp)) <= 0){
                        // query is bigger or equal  to center substring
                        l = c+1; //right bisect
                        lLcp = getLCP(getSubstring(prefixList.get(l), 0, queryLength), query);
                    } else{
                        //query is bigger than center
                        r = c;// right bisect
                        rLcp = getLCP(getSubstring(prefixList.get(r), 0, queryLength), query);
                    }
                }
                // convert prefixtable hits to actual genome indices
                leftBound = prefixTable.get(prefixList.get(leftBound)).getIndices(0);
                r = prefixTable.get(prefixList.get(r)).getIndices(1);
                for (int i = leftBound; i<r; i++){
                    hits.add(suffixArray.get(i).intValue());
                }   
                //stop timing
                long endTime = System.nanoTime();
                long duration = (endTime-startTime);
                System.out.println("Query of "+ query.length() + " length took " + duration);
                return hits;
                
            } else { 
                //query is longer than prefix, proceed with standard search
                String queryPrefix = query.substring(0, prefixLength);
                if (prefixTable.get(queryPrefix)==null){
                    long endTime = System.nanoTime();
                    long duration = (endTime-startTime);
                    System.out.println("Query of "+ query.length() + " length took " + duration);
                    return hits; //no matches if prefix can not be found
                } else{
                    l = prefixTable.get(queryPrefix).getIndices(0);
                    r = prefixTable.get(queryPrefix).getIndices(1);
                }             
            }
        }
        //binary search for lower bound
        int lLcp = getLCP(getSubstring(genome, suffixArray.get(l).intValue(),queryLength), query);
        int rLcp = getLCP(getSubstring(genome, suffixArray.get(r).intValue(),queryLength), query);
        while (l< r){
            int minLcp = Math.min(lLcp, rLcp);
            c = (int)Math.floor((l+r)/2);
            //check middle shifting by minLCP number of characters
            if (getSubstring(genome, suffixArray.get(c).intValue()+minLcp,queryLength-minLcp).compareTo(query.substring(minLcp)) >= 0){
                // query is  smaller or equal  to center substring
                r = c; //left bisect
                rLcp = getLCP(getSubstring(genome, suffixArray.get(r).intValue(),queryLength), query);
            } else{
                //query is bigger than center
                l = c+1;// right bisect
                lLcp = getLCP(getSubstring(genome, suffixArray.get(l).intValue(),queryLength), query);
            }
        }
        int leftBound = l;
        //binary search for upper bound (might as well start at the left bound)
        r = suffixArray.size()-1;
        rLcp = getLCP(getSubstring(genome, suffixArray.get(r).intValue(),queryLength), query);
        while (l<r){
            int minLcp = Math.min(lLcp, rLcp);
            c = (int)Math.floor((l+r)/2);
            //check middle
            if (getSubstring(genome, suffixArray.get(c).intValue()+minLcp,queryLength-minLcp).compareTo(query.substring(minLcp)) <= 0){
                // query is bigger or equal  to center substring
                l = c+1; //right bisect
                lLcp = getLCP(getSubstring(genome, suffixArray.get(l).intValue(),queryLength), query);
            } else{
                //query is bigger than center
                r = c;// right bisect
                rLcp = getLCP(getSubstring(genome, suffixArray.get(r).intValue(),queryLength), query);
            }
        }
        // convert suffix Array hits to actual genome indices
        for (int i = leftBound; i<r; i++){
            hits.add(suffixArray.get(i).intValue());
        }   
        //stop timing
        long endTime = System.nanoTime();
        long duration = (endTime-startTime);
        System.out.println("Query of "+ query.length() + " length took " + duration);
        return hits;
    }

    public static void main(String[] args) {
        String binPath, queryPath, queryMode, output;
        if (args.length==4){
            binPath = args[0]; //path to binary file
            queryPath = args[1]; //path to query fasta file
            queryMode = args[2]; //naive or simpaccel
            output = args[3]; //file path for output
        } else {
            System.err.println("Invalid arguments count: " + args.length);
            return;
        }
        // read in protobuf suffix array
        SuffixArrayProto.SuffixMessage suffixMessage = ReadSABin(binPath);
        String genome = suffixMessage.getGenome();
        System.out.println("Reference is " + genome.length() + " characters long");
        List<Integer> suffixArray = suffixMessage.getSuffixArrayList(); //suffixArray.get(i).intValue()
        Map<String,SuffixArrayProto.indexInterval> prefixTable = suffixMessage.getPreftabMap(); //prefixTable.get('pattern').getIndices(index)
        if (prefixTable.size()==0){
           System.out.println("No prefix table found nor used");
        } else {
            System.out.println("Prefix table found and used");
        }
       
        // read in fasta query
        Map<String,String> queries = ReadFastaQuery(queryPath);
        //conduct query and write output to tab seperated file
        try(FileWriter outputFile = new FileWriter(output)){
            if (queryMode.equals("naive")){
                System.out.println("Using naive query");
                for (Map.Entry<String, String> entry : queries.entrySet()){
                    List<Integer> queryResults = naiveQuery(entry.getValue(), genome, suffixArray, prefixTable);
                    String outputLine = entry.getKey() + "\t" + queryResults.size() + "\t" + 
                        queryResults.stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n";
                    outputFile.write(outputLine);
                }
                outputFile.close();
            }else if (queryMode.equals("simpaccel")){
                System.out.println("Using simple accelerant query");
                for (Map.Entry<String, String> entry : queries.entrySet()){
                    List<Integer> queryResults = naiveQuery(entry.getValue(), genome, suffixArray, prefixTable);
                    String outputLine = entry.getKey() + "\t" + queryResults.size() + "\t" + 
                        queryResults.stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n";
                    outputFile.write(outputLine);
                }
                outputFile.close();
            } else{
                outputFile.close();
                System.out.println("invalid querymode presented, expecting naive or simpaccel");
                System.exit(1);  
            }
        } catch (IOException e){
            System.out.println("Files failed to write");
            System.exit(1);
        }
        
        System.out.println("All queries complete and saved");
       
    }
}
