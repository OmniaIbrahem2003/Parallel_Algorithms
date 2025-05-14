#define _CRT_SECURE_NO_WARNINGS
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <algorithm>
#include <ctime>

using namespace std;

// Comparison function for qsort (used in Sample Sort)
static int intcompare(const void* i, const void* j) {
    int a = *(const int*)i;
    int b = *(const int*)j;
    return (a > b) - (a < b);
}

// Prime Number Finder
bool is_prime(long long n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (long long i = 5; i * i <= n; i += 6)
        if (n % i == 0 || n % (i + 2) == 0)
            return false;
    return true;
}

void runPrimeFinder(int rank, int size) {
    vector<long long> data;
    string input_file, output_file;
    double t_start, t_end;

    if (rank == 0) {
        cout << "=====================================\n";
        cout << "  PARALLEL PRIME NUMBER FINDER\n";
        cout << "=====================================\n";
        cout << "Processes: " << size << "\n\n";

        cout << "Enter input file path: ";
        cin >> input_file;
        cout << "Enter output file path: ";
        cin >> output_file;

        ifstream in(input_file);
        if (!in) {
            cerr << "ERROR: Cannot open input file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        string line;
        while (getline(in, line)) {
            stringstream ss(line);
            long long num;
            while (ss >> num) {
                data.push_back(num);
            }
        }
        in.close();

        if (data.empty()) {
            cerr << "ERROR: No numbers found in input file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        cout << "\nNumbers to check: " << data.size() << "\n";
        t_start = MPI_Wtime();
    }

    int total_count;
    if (rank == 0) total_count = data.size();
    MPI_Bcast(&total_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int chunk = total_count / size;
    vector<long long> subset(chunk);

    MPI_Scatter(data.data(), chunk, MPI_LONG_LONG, subset.data(), chunk, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        int remain = total_count % size;
        subset.insert(subset.end(), data.begin() + chunk * size, data.end());
    }

    unordered_set<long long> primes_local;
    for (long long n : subset) {
        if (is_prime(n)) {
            primes_local.insert(n);
        }
    }

    vector<long long> primes_vec(primes_local.begin(), primes_local.end());

    if (rank == 0) {
        unordered_set<long long> primes_global;
        primes_global.insert(primes_vec.begin(), primes_vec.end());

        for (int src = 1; src < size; src++) {
            int count;
            MPI_Recv(&count, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            vector<long long> recv_buf(count);
            MPI_Recv(recv_buf.data(), count, MPI_LONG_LONG, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            primes_global.insert(recv_buf.begin(), recv_buf.end());
        }

        t_end = MPI_Wtime();
        double elapsed = t_end - t_start;

        vector<long long> sorted_primes(primes_global.begin(), primes_global.end());
        sort(sorted_primes.begin(), sorted_primes.end());

        ofstream out(output_file);
        if (!out) {
            cerr << "ERROR: Cannot open output file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        out << "Unique prime numbers found (" << sorted_primes.size() << "):\n";
        for (long long p : sorted_primes) {
            out << p << "\n";
        }
        out.close();

        cout << "\nRESULTS\n";
        cout << "-------\n";
        cout << "Unique primes found: " << sorted_primes.size() << "\n";
        cout << fixed << setprecision(4);
        cout << "Execution time: " << elapsed << " seconds\n";
        cout << "Results saved to: " << output_file << "\n";
        cout << "---------------------------------\n";
    }
    else {
        int count = primes_vec.size();
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(primes_vec.data(), count, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    }
}

// Parallel Search Target
void runParallelSearch(int rank, int size) {
    vector<int> all_numbers;
    int target;
    string inputPath, outputPath;
    double start_time, end_time;

    if (rank == 0) {
        cout << "=====================================\n";
        cout << "        PARALLEL SEARCH TARGET       \n";
        cout << "=====================================\n";
        cout << "Processes: " << size << "\n\n";

        cout << "Enter input file path: ";
        cin >> inputPath;
        cout << "Enter target number: ";
        cin >> target;
        cout << "Enter output file path: ";
        cin >> outputPath;

        ifstream infile(inputPath);
        if (!infile) {
            cerr << "ERROR: Cannot open input file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int num;
        while (infile >> num) {
            all_numbers.push_back(num);
        }
        infile.close();

        if (all_numbers.empty()) {
            cerr << "ERROR: No numbers found in input file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        cout << "\nNumbers to search: " << all_numbers.size() << "\n";
        start_time = MPI_Wtime();
    }

    int total_numbers = all_numbers.size();
    MPI_Bcast(&total_numbers, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&target, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        all_numbers.resize(total_numbers);
    }
    MPI_Bcast(all_numbers.data(), total_numbers, MPI_INT, 0, MPI_COMM_WORLD);

    int chunkSize = total_numbers / size;
    int startIdx = rank * chunkSize;
    int endIdx = (rank == size - 1) ? total_numbers : startIdx + chunkSize;

    int foundIndex = -1;
    for (int i = startIdx; i < endIdx; ++i) {
        if (all_numbers[i] == target) {
            foundIndex = i;
            break;
        }
    }

    vector<int> foundIndices(size);
    MPI_Gather(&foundIndex, 1, MPI_INT, foundIndices.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    int globalIndex = -1;
    MPI_Reduce(&foundIndex, &globalIndex, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        end_time = MPI_Wtime();

        ofstream outfile(outputPath, ios::app);
        if (!outfile) {
            cerr << "ERROR: Cannot open output file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        time_t now = time(0);
        char* dt = ctime(&now);

        outfile << "===========================================\n";
        outfile << "Run Timestamp: " << dt;
        outfile << "Target: " << target << "\n";
        outfile << "Number of Processes: " << size << "\n\n";

        for (int i = 0; i < size; ++i) {
            if (foundIndices[i] != -1)
                outfile << "Process " << i << ": Found at global index " << foundIndices[i] << "\n";
            else
                outfile << "Process " << i << ": Not found in its chunk.\n";
        }

        if (globalIndex != -1) {
            int finderRank = -1;
            for (int i = 0; i < size; ++i) {
                if (foundIndices[i] == globalIndex) {
                    finderRank = i;
                    break;
                }
            }
            outfile << "\nResult: Value found at global index " << globalIndex
                << " by process " << finderRank << "\n";
        }
        else {
            outfile << "\nResult: Value not found in array.\n";
        }

        outfile << fixed << setprecision(6);
        outfile << "Execution Time: " << (end_time - start_time) << " seconds\n";
        outfile << "-------------------------------------------\n\n";

        outfile.close();

        cout << "\nResult written to: " << outputPath << "\n";
        cout << fixed << setprecision(6);
        cout << "Execution Time: " << (end_time - start_time) << " seconds\n";
    }
}

// Radix Sort
int findMax(const vector<int>& values) {
    return *max_element(values.begin(), values.end());
}

void stableDigitSort(vector<int>& segment, int exponent) {
    int n = segment.size();
    vector<int> sorted(n);
    int count[10] = { 0 };

    for (int i = 0; i < n; i++)
        count[(segment[i] / exponent) % 10]++;

    for (int i = 1; i < 10; i++)
        count[i] += count[i - 1];

    for (int i = n - 1; i >= 0; i--) {
        int digit = (segment[i] / exponent) % 10;
        sorted[count[digit] - 1] = segment[i];
        count[digit]--;
    }

    for (int i = 0; i < n; i++)
        segment[i] = sorted[i];
}

void radixSort(vector<int>& subarray, int maximum) {
    for (int exponent = 1; maximum / exponent > 0; exponent *= 10)
        stableDigitSort(subarray, exponent);
}

void runRadixSort(int rank, int size) {
    vector<int> inputNumbers;
    string inputFile, outputFile;
    double timerStart, timerEnd;

    if (rank == 0) {
        cout << "=====================================\n";
        cout << "         RADIX SORT      \n";
        cout << "=====================================\n";
        cout << "Total MPI Processes: " << size << "\n\n";

        cout << "Enter path to input file: ";
        cin >> inputFile;
        cout << "Enter path to output file: ";
        cin >> outputFile;

        ifstream fin(inputFile);
        if (!fin) {
            cerr << "ERROR: Unable to open input file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int value;
        while (fin >> value)
            inputNumbers.push_back(value);
        fin.close();

        if (inputNumbers.empty()) {
            cerr << "ERROR: No data found in the input file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (inputNumbers.size() % size != 0) {
            cerr << "ERROR: Input size must be divisible by process count.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        cout << "\nInput size: " << inputNumbers.size() << "\n";
        timerStart = MPI_Wtime();
    }

    int totalSize;
    if (rank == 0) totalSize = inputNumbers.size();
    MPI_Bcast(&totalSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int segmentSize = totalSize / size;
    vector<int> localChunk(segmentSize);

    MPI_Scatter(inputNumbers.data(), segmentSize, MPI_INT,
        localChunk.data(), segmentSize, MPI_INT,
        0, MPI_COMM_WORLD);

    int localPeak = findMax(localChunk);
    int absoluteMax;
    MPI_Allreduce(&localPeak, &absoluteMax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    radixSort(localChunk, absoluteMax);

    if (rank == 0) inputNumbers.clear();
    inputNumbers.resize(totalSize);

    MPI_Gather(localChunk.data(), segmentSize, MPI_INT,
        inputNumbers.data(), segmentSize, MPI_INT,
        0, MPI_COMM_WORLD);

    if (rank == 0) {
        timerEnd = MPI_Wtime();
        ofstream fout(outputFile);
        if (!fout) {
            cerr << "ERROR: Unable to open output file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (size_t i = 0; i < inputNumbers.size(); ++i) {
            fout << inputNumbers[i];
            if (i != inputNumbers.size() - 1)
                fout << " ";
        }
        fout << "\n";
        fout.close();

        cout << "\n=== SORT COMPLETE ===\n";
        cout << "Output written to: " << outputFile << "\n";
        cout << fixed << setprecision(4);
        cout << "Elapsed time: " << (timerEnd - timerStart) << " seconds\n";
        cout << "=============================\n";
    }
}

// Bitonic Sort
bool isPowerOfTwo(int number) {
    return (number > 0) && ((number & (number - 1)) == 0);
}

void performBitonicMerge(vector<int>& localSegment, int segmentSize, int rank, int partnerRank, bool ascending, MPI_Comm comm) {
    vector<int> receivedSegment(segmentSize);
    MPI_Sendrecv(localSegment.data(), segmentSize, MPI_INT, partnerRank, 0,
        receivedSegment.data(), segmentSize, MPI_INT, partnerRank, 0,
        comm, MPI_STATUS_IGNORE);

    vector<int> mergedSegment(2 * segmentSize);
    copy(localSegment.begin(), localSegment.end(), mergedSegment.begin());
    copy(receivedSegment.begin(), receivedSegment.end(), mergedSegment.begin() + segmentSize);

    sort(mergedSegment.begin(), mergedSegment.end());

    if ((ascending && rank < partnerRank) || (!ascending && rank > partnerRank)) {
        copy(mergedSegment.begin(), mergedSegment.begin() + segmentSize, localSegment.begin());
    }
    else {
        copy(mergedSegment.begin() + segmentSize, mergedSegment.end(), localSegment.begin());
    }
}

void runBitonicSort(int rank, int size) {
    if (rank == 0 && !isPowerOfTwo(size)) {
        cerr << "ERROR: The number of processes must be a power of 2 for Bitonic Sort." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    vector<int> allData;
    int totalElements = 0;
    string inputFile, outputFile;
    double startTime, endTime;

    if (rank == 0) {
        cout << "===============================================\n";
        cout << "                    Bitonic Sort               \n";
        cout << "===============================================\n";

        cout << "Enter the path to the input file: ";
        cin >> inputFile;
        cout << "Enter the path for the output file: ";
        cin >> outputFile;

        ifstream countStream(inputFile);
        if (!countStream.is_open()) {
            cerr << "ERROR: Cannot open input file." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int tempVal;
        while (countStream >> tempVal)
            totalElements++;
        countStream.close();

        cout << "\nTotal number of elements in input file: " << totalElements << endl;

        if (totalElements % size != 0) {
            cerr << "ERROR: Number of elements must be divisible by number of processes." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        allData.resize(totalElements);
        ifstream dataStream(inputFile);
        for (int i = 0; i < totalElements; ++i)
            dataStream >> allData[i];
        dataStream.close();
    }

    MPI_Bcast(&totalElements, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int localSize = totalElements / size;
    vector<int> localData(localSize);

    MPI_Scatter(allData.data(), localSize, MPI_INT,
        localData.data(), localSize, MPI_INT,
        0, MPI_COMM_WORLD);

    startTime = MPI_Wtime();

    sort(localData.begin(), localData.end());

    for (int k = 2; k <= size; k *= 2) {
        for (int j = k / 2; j > 0; j /= 2) {
            int partnerRank = rank ^ j;
            bool sortAscending = ((rank & k) == 0);
            performBitonicMerge(localData, localSize, rank, partnerRank, sortAscending, MPI_COMM_WORLD);
        }
    }

    if (rank == 0) allData.resize(totalElements);
    MPI_Gather(localData.data(), localSize, MPI_INT,
        allData.data(), localSize, MPI_INT,
        0, MPI_COMM_WORLD);

    if (rank == 0) {
        endTime = MPI_Wtime();

        ofstream outStream(outputFile);
        if (!outStream.is_open()) {
            cerr << "ERROR: Cannot open output file." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < totalElements; ++i)
            outStream << allData[i] << " ";
        outStream << endl;
        outStream.close();

        cout << "\nSorted output written to: " << outputFile << endl;
        cout << "Execution Time: " << (endTime - startTime) << " seconds\n";
    }
}

// Sample Sort 
void runBucketSort(int rank, int size) {
    int NoofElements = 0, NoofElements_Bloc, NoElementsToSort;
    int* Input = nullptr, * InputData = nullptr;
    int* Splitter = nullptr, * AllSplitter = nullptr;
    int* Buckets = nullptr, * BucketBuffer = nullptr, * LocalBucket = nullptr;
    int* OutputBuffer = nullptr, * Output = nullptr;
    double start_time, end_time;
    string inputFile, outputFile;
    const int Root = 0;

    if (rank == Root) {
        cout << "=====================================\n";
        cout << "         SAMPLE SORT       \n";
        cout << "=====================================\n";
        cout << "Processes: " << size << "\n\n";

        cout << "Enter the path to the input file: ";
        cin >> inputFile;
        cout << "Enter the path for the output file: ";
        cin >> outputFile;

        // Read input file to count elements
        ifstream countStream(inputFile);
        if (!countStream.is_open()) {
            cerr << "ERROR: Cannot open input file!\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int tempVal;
        vector<int> inputNumbers;
        while (countStream >> tempVal)
            inputNumbers.push_back(tempVal);
        countStream.close();

        NoofElements = inputNumbers.size();
        if (NoofElements < 3) {
            cerr << "ERROR: Number of elements must be at least 3\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (NoofElements % size != 0) {
            cerr << "ERROR: Number of elements must be divisible by number of processes\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        Input = new int[NoofElements];
        Output = new int[NoofElements];
        if (Input == nullptr || Output == nullptr) {
            cerr << "Root process: Memory allocation failed for Input or Output\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < NoofElements; ++i) {
            Input[i] = inputNumbers[i];
        }
    }

    start_time = MPI_Wtime();

    // Broadcast the number of elements
    MPI_Bcast(&NoofElements, 1, MPI_INT, Root, MPI_COMM_WORLD);

    NoofElements_Bloc = NoofElements / size;
    InputData = new int[NoofElements_Bloc];
    if (InputData == nullptr) {
        cerr << "Process " << rank << ": Memory allocation failed for InputData\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Scatter(Input, NoofElements_Bloc, MPI_INT, InputData, NoofElements_Bloc, MPI_INT, Root, MPI_COMM_WORLD);

    // Sort local data
    qsort(InputData, NoofElements_Bloc, sizeof(int), intcompare);

    // Gather locally sorted data (optional, used later)
    MPI_Gather(InputData, NoofElements_Bloc, MPI_INT, Output, NoofElements_Bloc, MPI_INT, Root, MPI_COMM_WORLD);

    // Choose local splitters
    Splitter = new int[size - 1];
    if (Splitter == nullptr) {
        cerr << "Process " << rank << ": Memory allocation failed for Splitter\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < size - 1; i++) {
        Splitter[i] = InputData[(NoofElements_Bloc / size) * (i + 1)];
    }

    // Gather all local splitters at root
    if (rank == Root) {
        AllSplitter = new int[size * (size - 1)];
        if (AllSplitter == nullptr) {
            cerr << "Root process: Memory allocation failed for AllSplitter\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Gather(Splitter, size - 1, MPI_INT, AllSplitter, size - 1, MPI_INT, Root, MPI_COMM_WORLD);

    // Root sorts all splitters and selects global splitters
    if (rank == Root) {
        qsort(AllSplitter, size * (size - 1), sizeof(int), intcompare);
        for (int i = 0; i < size - 1; i++) {
            Splitter[i] = AllSplitter[(size - 1) * (i + 1)];
        }
    }

    // Broadcast global splitters
    MPI_Bcast(Splitter, size - 1, MPI_INT, Root, MPI_COMM_WORLD);

    // Create buckets
    Buckets = new int[(NoofElements_Bloc + 1) * size];
    if (Buckets == nullptr) {
        cerr << "Process " << rank << ": Memory allocation failed for Buckets\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < size; i++) {
        Buckets[i * (NoofElements_Bloc + 1)] = 0;
    }

    int j = 0, k = 1;
    for (int i = 0; i < NoofElements_Bloc; i++) {
        if (j < size - 1 && InputData[i] >= Splitter[j]) {
            Buckets[j * (NoofElements_Bloc + 1)] = k - 1;
            j++;
            k = 1;
            i--;
        }
        else {
            Buckets[j * (NoofElements_Bloc + 1) + k] = InputData[i];
            k++;
        }
    }
    Buckets[j * (NoofElements_Bloc + 1)] = k - 1;

    // Allocate bucket buffer
    BucketBuffer = new int[(NoofElements_Bloc + 1) * size];
    if (BucketBuffer == nullptr) {
        cerr << "Process " << rank << ": Memory allocation failed for BucketBuffer\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // All-to-all communication
    MPI_Alltoall(Buckets, NoofElements_Bloc + 1, MPI_INT, BucketBuffer, NoofElements_Bloc + 1, MPI_INT, MPI_COMM_WORLD);

    // Rearrange received data
    LocalBucket = new int[2 * NoofElements_Bloc];
    if (LocalBucket == nullptr) {
        cerr << "Process " << rank << ": Memory allocation failed for LocalBucket\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int count = 1;
    for (j = 0; j < size; j++) {
        int c = BucketBuffer[j * (NoofElements_Bloc + 1)];
        for (k = 1; k <= c; k++) {
            LocalBucket[count++] = BucketBuffer[j * (NoofElements_Bloc + 1) + k];
        }
    }
    LocalBucket[0] = count - 1;

    // Sort local bucket
    NoElementsToSort = LocalBucket[0];
    qsort(&LocalBucket[1], NoElementsToSort, sizeof(int), intcompare);

    // Prepare buffers for gathering
    if (rank == Root) {
        OutputBuffer = new int[2 * NoofElements];
        Output = new int[NoofElements];
        if (OutputBuffer == nullptr || Output == nullptr) {
            cerr << "Root process: Memory allocation failed for Output buffers\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Gather(LocalBucket, 2 * NoofElements_Bloc, MPI_INT, OutputBuffer, 2 * NoofElements_Bloc, MPI_INT, Root, MPI_COMM_WORLD);

    if (rank == Root) {
        count = 0;
        for (j = 0; j < size; j++) {
            int c = OutputBuffer[j * 2 * NoofElements_Bloc];
            for (int i = 1; i <= c; i++) {
                Output[count++] = OutputBuffer[j * 2 * NoofElements_Bloc + i];
            }
        }
        end_time = MPI_Wtime();

        ofstream fp(outputFile);
        if (!fp.is_open()) {
            cerr << "Root process: Cannot open output file\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fp << "Number of Elements to be sorted: " << NoofElements << "\n";
        fp << "The sorted sequence is: ";
        for (int i = 0; i < NoofElements; i++) {
            fp << Output[i];
            if (i < NoofElements - 1) fp << " ";
        }
        fp << "\n";
        fp << "Number of Elements to be sorted: " << NoofElements << "\n";
        fp << "Number of processors working together: " << size << "\n";
        fp << "Elapsed Time: " << (end_time - start_time) << " seconds\n";
        fp.close();

        cout << "Elapsed Time: " << (end_time - start_time) << " seconds\n";
    }

    // Free memory
    delete[] InputData;
    delete[] Splitter;
    delete[] Buckets;
    delete[] BucketBuffer;
    delete[] LocalBucket;
    if (rank == Root) {
        delete[] Input;
        delete[] OutputBuffer;
        delete[] Output;
        delete[] AllSplitter;
    }
}

// Main function with single list-based prompt
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int choice = 0;
    if (rank == 0) {
        cout << "=====================================\n";
        cout << "  MPI PARALLEL PROGRAMS SELECTION\n";
        cout << "=====================================\n";
        cout << "Number of Processes: " << size << "\n\n";
        cout << "Select a program to run:\n";
        cout << "1. Prime Number Finder\n";
        cout << "2. Parallel Search Target\n";
        cout << "3. Radix Sort\n";
        cout << "4. Bitonic Sort\n";
        cout << "5. Sample Sort\n";
        cout << "Enter your choice (1-5): ";
        cin >> choice;
    }

    // Broadcast the choice
    MPI_Bcast(&choice, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Execute the selected program
    switch (choice) {
    case 1:
        runPrimeFinder(rank, size);
        break;
    case 2:
        runParallelSearch(rank, size);
        break;
    case 3:
        runRadixSort(rank, size);
        break;
    case 4:
        runBitonicSort(rank, size);
        break;
    case 5:
        runBucketSort(rank, size);
        break;
    default:
        if (rank == 0) {
            cerr << "ERROR: Invalid choice! Please select a number between 1 and 5.\n";
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
        cout << "\nProgram completed.\n";
    }

    MPI_Finalize();
    return 0;
}