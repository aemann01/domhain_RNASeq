using System;
using System.IO;
using System.Linq;
using System.Collections.Concurrent;
using System.Threading.Tasks;
using System.Text.RegularExpressions;
using System.Threading;

class Program
{
    static void Main()
    {
        var geneIds = File.ReadLines("clean_annots.txt")
            .Select(line => line.Split()[1].Replace("tag/", "Geneid/"))
            .ToList();

        var readLines = File.ReadAllLines("read_counts2.txt");
        var output = new ConcurrentBag<string>();

        // Track progress
        int processedCount = 0;
        int totalGeneIds = geneIds.Count;

        // Parallel processing with progress tracking
        Parallel.ForEach(geneIds, new ParallelOptions { MaxDegreeOfParallelism = 190 }, geneId =>
        {
            var pattern = new Regex($@"\b{Regex.Escape(geneId)}\b");
            foreach (var line in readLines)
            {
                if (pattern.IsMatch(line))
                {
                    output.Add(line);
                    break;
                }
            }

            // Thread-safe progress update
            int currentProcessed = Interlocked.Increment(ref processedCount);

            // Display progress every 1000 gene IDs processed
            if (currentProcessed % 1000 == 0 || currentProcessed == totalGeneIds)
            {
                Console.WriteLine($"Progress: {currentProcessed} of {totalGeneIds} gene IDs processed.");
            }
        });

        // Write output to file
        File.WriteAllLines("red_counts.txt", output.OrderBy(x => x)); // Ordering is optional
        Console.WriteLine("Processing complete. Results written to gene_counts.txt");
    }
}
