//
// This file contains code for editing and merging some intermediate files.
// This file is not used in diff matrix calculations.
//

#include <direct.h>
#include <assert.h>
#include "Set.h"
#include "Constants.h"
#include "Globals.h"

void MergeMawRC()
{
    int i;
    Set original, rc, merged;
    char strFileName[50];
    
    for (i = 0; i < g_numGenes; i++)
    {
        sprintf(strFileName, "NoRC\\%s.maw.txt", g_strSpeciesFullName[i]);
        original = Set::CreateFromFile(strFileName);

        sprintf(strFileName, "WithRC\\%s.maw.txt", g_strSpeciesFullName[i]);
        rc = Set::CreateFromFile(strFileName);

        merged = original.Union(rc);

        //
        // Write the merged set out to file
        //
        sprintf(strFileName, "%s.raw.txt", g_strSpeciesFullName[i]);
        freopen(strFileName, "w", stdout); 
        printf(">1\n");
        merged.Print();
        fclose(stdout);
    }
}

void GetEditedRawFilePath(int i, int j, char strPathOut[])
{
    sprintf(
        strPathOut,
        "%s\\%s_%s.fasta",
        g_strSpeciesFullName[i], 
        g_strSpeciesShortName[i], 
        g_strSpeciesFullName[j]
        );

    mkdir(g_strSpeciesFullName[i]);
}

void GetRawFilePath(int i, int j, char *strBaseDir, char strPathOut[])
{
    char strFile[50];

    sprintf(
        strFile,
        "%s\\%s_%s.raw.txt",
        g_strSpeciesFullName[i], 
        g_strSpeciesShortName[i], 
        g_strSpeciesFullName[j]
        );

    sprintf(strPathOut, "%s\\%s", strBaseDir, strFile);
}

void EditRaw()
{
    //
    // Raw files contained a leading decimal number in each line
    // Removing that
    //
    int i, j;
    char strSourcePath[50], strDestPath[50];
    char buf[MAX_WORD_LEN];
    
    for (i = 0; i < g_numGenes; i++)
    {
        for (j = 0; j < g_numGenes; j++)
        {
            if (i == j)
            {
                continue;
            }

            //
            // Make the file name from the gene indices
            //
            GetRawFilePath(i, j, "Original", strSourcePath);
            GetEditedRawFilePath(i, j, strDestPath);
            
            freopen(strSourcePath, "r", stdin);
            freopen(strDestPath, "w", stdout); 

            while (EOF != scanf("%*d%s", buf))
            {
                printf("%s\n", buf);
            }
            printf("\n");
            
            fclose(stdin);
            fclose(stdout);
        }
    }
}

void MergeRawRC()
{
    int i, j;
    Set original, rc, merged;
    char strFullPath[50];
    
    for (i = 0; i < g_numGenes; i++)
    {
        for (j = 0; j < g_numGenes; j++)
        {
            if (i == j)
            {
                continue;
            }

            //
            // Make the file name from the gene indices
            //
            GetRawFilePath(i, j, "NoRC", strFullPath);
            original = Set::CreateFromFile(strFullPath);

            GetRawFilePath(i, j, "WithRC", strFullPath);
            rc = Set::CreateFromFile(strFullPath);

            merged = original.Union(rc);

            //
            // Write the merged set out to file
            //
            GetEditedRawFilePath(i, j, strFullPath);
            freopen(strFullPath, "w", stdout); 
            merged.Print();
            fclose(stdout);
        }
    }
}