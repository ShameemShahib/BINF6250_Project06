# Introduction
In this project, our group implements the Neighbor-Joining (NJ) algorithm for phylogenetic tree construction using Smith-Waterman alignment scores. Neighbor-Joining builds trees by iteratively joining pairs of sequences based on a distance matrix, producing an unrooted tree that reflects evolutionary relationships. To generate these distances, we use the Smith-Waterman algorithm, which computes optimal local alignments between sequences and captures biologically meaningful similarities. By converting these alignment scores into distances and applying the NJ algorithm, we can reconstruct relationships between sequences from real biological data. This approach is especially useful for analyzing datasets with varying evolutionary rates and complex sequence similarities.

# Pseudocode
Put pseudocode in this box:

```
Read_fasta function:
Create an empty dict for sequences
Create an empty string for our current_header and an empty string for our current_sequence
Open the file
Iterate through the lines of the file
Check if the line is header (“>” at start) or a sequence
	if it’s a header
		We need to pack up the string to which we’ve been appending lines from the current sequence and add it to the dictionary 
			This requires us to “remember” what header that sequence was associated with (don’t update our current_header variable until after packing up this sequence)
			Here, we also need to check if we even have a header or a sequence that we’ve already been working with, as this will run on the first header of the file as well. So only add the entry if we have a header and a sequence
		Then we update what our current_header is and clear our current_sequence
	
if the current line isn’t a header, we just append the line to current_sequence variable (strip it first to remove newlines)

We just return the dictionary at the end

Smith_waterman function: 
import the smith_waterman function from the textdistance library :D (it constructs the alignment matrix for us and calculates the similarity or distance score)
	build the alignment matrix for the two sequences using the scoring params we’re given here
	Grab the highest value in the alignment matrix → this is our similarity score, as that score’s position is where we’d start the traceback from to complete the local alignment
Return the similarity score

Build_distance_matrix function: 
Initialize an N x N matrix, where N is number of sequences, and the value in each position is 0
For every combination of two sequences (nested for loop)
	Run smith_waterman to get the distance score
	Find the position representing the two sequences we’re looking at (this is just based on their index in the list, and we can iterate thru the list with enumerate so we have this ready to go), and update the value to the distance score we just calculated

Return the distance matrix


Graph object: 


Node object:
Store the seqID
Store branch length
Store the node’s parent

neighbor_joining function: 
make a working copy of the distance matrix
initialize our graph object
pick 3 nodes:
	Get the two that are the largest distance apart from each other (combo of seq’s on the distance matrix with the highest value), these are our “candidates”
	pick a random third sequence from the remaining on the matrix
	Determine which of the two candidates are less related to the third sequence, and that sequence is our sequence that we will bald and trim
Make a node object for this sequence:
calculate the branch length for the sequence that we’re going to bald
subtract the branch length from all the distances for the sequence in question, copy these values into a vector of balded distances for this sequence
Remove the row and column associated with this sequence

Call neighbor_joining if the pruned matrix is larger than 2x2
If the matrix is 2x2, we make a node object for each of the two remaining sequences:
	Their branch lengths will both be equal to their distance to one another
	Add just these two nodes to the graph

Now that we’re breaking out of the recursion, we need to: 
make an internal node between the last two nodes we touched
NOTE: the most recent leaf (the one that got trimmed off during this recursive layer) will SHARE a parent with the internal node we’re making
	if we’re just starting the graph:
The last two nodes we touched will just be the two leaves remaining on the distance matrix
	If we already have nodes on the graph:
One of these nodes is going to be the leaf that we trimmed off in the previous function call/recursive layer/step/whatever
The other one would be the last internal node that got made, in our newick representation would just be sort of a list of all the nodes downstream of it
Now we have the distance between the parent of the current leaf to each of the last two nodes touched (the previous leaf’s parent and the previous internal node)
Calculate the branch length for the previous internal node USING the parent of the leaf that got trimmed off in this function call/recursive layer, this tells us the distance between these two and the new internal node we’re making right now
NOTE: we use the balded distances as a representation of how distant a node’s parent is from every other leaf, so we can determine how distant the new leaf’s parent is from the internal node we’re making

So now we add that internal node to the graph + the current leaf (this is actually the last thing the function would need to do, since we won’t need to make an internal node between the last leaf to be trimmed and the final internal node (the one that acts like a root to the rest of the leaves on the graph) since our graph is unrooted)


Example, we have the last two leaves in the matrix: A and B
distance from A to B is the branch length for both A and B
We introduce leaf C
Leaf C has a branch length, and that is the distance between C and Parent-C
Now A and B also share a parent: AB
Parent-C is also the parent to AB and C, so we can call it ABC
Now we know the distance from:
	C to ABC
	A to B (which is A to AB + B to AB)
	ABC to A
	ABC to B
So we use these distances: A→B, ABC→A, ABC→B to determine the following distances:
A to AB
B to AB
ABC to AB
Now we can update the branch lengths for A and B, so that A’s branch length is A→AB and B’s is B→AB
Additionally, the branch length for AB is equal to the distance from ABC to AB
but ABC doesn’t have a branch length

Now we add a new leaf, D
D has a branch length that’s already been calculated, and that will be the distance between D and parent-D
Now we’ve already established that AB and C share the parent ABC
ABC and D share a parent, so parent-D can be called ABCD
Now we know the distance between parent-D (aka ABCD) and any other leaf (A, B, or C) because those are the balded distances we stored for node D
We also know how far D is from ABCD because that’s just D’s branch length
We’re looking for the distance from ABC to ABCD, which we can express using:
ABCD to C
|--> ABCD→C = ABCD→ABC + ABC→C
ABCD to either of A or B (let’s say it’s A)
|--> ABCD→A = ABCD→ABC + ABC→AB + AB→A
A to C
|--> A → C = C→ABC + ABC→AB + AB→A
So we rearrange algebraically and this is just our equation for branch length, but the reason I’m writing it out is because if we picked A and B instead of A and C, we would get:
(ABCD to B) + (ABCD to A) - (A to B)
(ABCD→ABC + ABC→AB + AB→A) + (ABCD→ABC + ABC→AB + AB→B) - (B→AB + AB→A)
(ABCD→ABC + ABC→AB) + (ABCD→ABC + ABC→AB)
2(ABC→AB) + 2(ABCD→ABC), which isn’t what we want
So now we have the distance from ABC to ABCD, this is the branch length for ABC\
ABCD’s branch length can’t be calculated yet until we introduce leaf E

I think this makes a rooted tree though, with the root being ABCD

Here’s how we make it unrooted:

Two nodes share a parent if the distance between them is equal to the sum of their branch lengths (pretty sure this is just always true)
e.g. D and C share a parent if [D→ parent-D] + [C→parent-C] is equal to [D→C]. If they share a parent, then D→C, which is normally [D→parent-D] + [parent-D→parent→C] + [parent-C→C], will have parent-D→parent-C = 0.  
So whenever we have a new leaf to place on the graph, we just make this assertion, that the new leaf shares a parent with whatever node we’re looking at, for every node on the graph (doesn’t matter if internal or terminal)
So we assert that the leaf shares a parent with some node on the graph, if they do, we have to check if the parent is already on the graph (e.g. if we add a leaf E to the example above that shares a parent with D). If that parent’s already on the graph, we just throw node E onto the graph at [E’s branch length] away from whatever node on the graph we decided was its parent
If it doesn’t exist, I think we do all that stuff above that I outlined with ABCD
```

# Successes
One of our biggest successes was our collaboration. We had effective communication through calls and texting, which made a huge difference when working through parts of the assignment. Being able to talk things out in real time helped us debug faster, share ideas, and stay on the same page. We found the textdistance library to be a valuable aid to our workflow, and we were able to help each other clear up doubts about the concepts involved in graph building. Overall, the preplanning phase helped build a solid foundation. Even when we embarked on work asynchronously, we were on the same page about what ideas we wanted to implement because of the prior planning.

# Struggles
As helpful as the pseudocode and planning was, the final neighbor joining function was also the most challenging part of this project. Before even coding, we had to spend time really understanding how the algorithm works and how to structure it properly, and we all had to spend a lot of time making sure we had the concept right so that we could put it into code correctly as well. Translating that plan into code was still tricky, especially when it came to setting up our objects. Implementing the textdistance library also brought some obstacles that we had to overcome, especially with how our graph was visualized. Overall, this project pushed us to be more methodical with both our planning and implementation, and showed how important it is to fully understand an algorithm before trying to code it.

## Reflection after initial incomplete

Our initial implementation of the neighbor-joining algorithm took the wrong approach to identifying which sequences to group at each step. Rather than finding the two most closely related nodes and merging them, we had been doing the inverse. We identified the node least related to the rest of the dataset and removed it from the matrix. While this sounds intuitively similar, it produces a fundamentally different tree and breaks the core logic of the algorithm.
The key insight that drove our revision was the Q-matrix. The Q-matrix is a transformation of the raw distance matrix that answers how much more related sequences i and j are to each other than each is to everything else. By using the Q-matrix to identify pairs, we simplified the process of determining which two nodes to merge at each step. The pair with the smallest Q-matrix value is the correct choice, and that clarity helped us revise our code to produce a tree with meaningful results
Beyond the algorithmic change, we also discovered several implementation bugs that had gone undetected. The tree structure was missing the ability to enumerate its own leaves, which likely contributed to difficulties printing the tree at the end. Additionally, we had an indexing error where the new distance matrix was being built with np.zeros(size+1, size+1) where size was the length of the matrix before the new node was added, but the code then wrote the new node's distances to index size-1 instead of size, which was always the wrong slot. Because size hadn't been updated after the resize, the index math was off by one throughout, leaving the last column consistently incorrect.
Through our revision, we were able to make significant changes to our code that ultimately helped us produce something that looked like a tree that could be interpreted. Considering the difficulties we faced in our intiial pass at the project, it was very satisfying to come back after learning where we went wrong and correct everything. We now not only feel more satisfied with our implementation, but we also feel like we have a much richer understanding of the algorithm as a whole after revisiting the project. 



# Personal Reflections
## Group Leader
Group leader's reflection on the project
Shameem Shahib:
This project was one of the more enjoyable and rewarding challenges for me so far this semester. We started with the foundational functions for reading in FASTA files and generating pairwise similarity scores, and when we got to Smith-Waterman, we made the practical decision to use the textdistance library's built-in implementation. The neighbor joining algorithm was definitely the most demanding piece, since translating the conceptual steps into working recursive code required a lot of careful thinking about how the matrix shrinks at each step and how to preserve enough information to build the tree back up afterwards. Having detailed notes to reference made the whole process a lot less intimidating and gave us a solid roadmap to follow when we got stuck. Our collaboration was a real strength throughout, as our meetings were consistently productive and we kept each other in the loop on progress, which meant no one was ever left behind or working in the dark. Overall, I'm really proud of how we came together as a group and the contributions we were all able to make.


## Other members
Aaronie Jersha Jenyfred: 
This was a challenging project in terms of understanding the complexity of algorithm and how each step contributes to building the tree. Translating the pseudocode for neighbour joining function into code was tricky especially managing the recursion part. While adding the leaf back, understanding how to use the last two nodes and compute the internal node’s branch length took some debugging. Additionally, understanding how graph structures were used to represent relationships between nodes was the toughest part for me.

Justin Wildman:
This was a very difficult project for me due to the time crunch from our group not being free to meet until Sunday and due how difficult it was to figure out how to decide which nodes are related to one another. The `textdistances` library was very convenient, but the time saved from not having to worry about doing the calculations for our distance matrix was more than eaten up by the time it took to implement the `Node` and `Tree_Graph` objects. It took quite a while to figure out what attributes the nodes and graph needed to have, whether a method would be more appropriate on one, the other, or called externally, and how to ensure that the branch lengths are "correct" (in the sense of not being an overestimate due to choosing two points of reference who are too closely related). Making the tree object's class methods be resilient against working with internal nodes and leaves was also much more annoying than originally anticipated, as getting the distances between an internal node and anything else is more involved than for leaves. I wish we had started this project earlier because it's a really cool algorithm and it would've been super interesting to see this through more completely. I think in my implementation of the node and tree classes, I blurred the line too much between what the `neighbor_joining()` function should be doing and what the objects need to handle, which I think led to a lot of redundancy in my code where it was already handled in `neighbor_joining()`. 

# Generative AI Appendix
As per the syllabus, none was used
