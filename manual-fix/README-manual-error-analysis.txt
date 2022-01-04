## README for manual error analysis and fixes

After the semi-automated effort to process the Fauci-email PDF, we performed an additional quality check and subsequent manual fix to correct a few errors in the automated processing we noticed that were due usually to challenges parsing through redacted information in the the sender/receiver/cc/subject information. This is a summary of the errors and how we fixed them in order to produced the released dataset.

Example errors
---------------

In a number of places, the cc information was pushed to the subject line, and the subject line was pushed to the body text.

In other places, the sender/receiver/cc information was correct, but the subject line was pushed to body, resulting in an empty subject line in our processed JSON file.


Overview of correction process
------------------------------

We corrected these types of errors in two batches:

* Errors in longer subject lines that looked suspicious (redacted, or appearing to have cc information in the subject line)
* Errors in emails that had an empty subject line (which could be indicative that the subject was processed as part of the email body text, and possibly other errors occurred as well in the sender/receiver/cc information)

We addressed these in steps, using the following code steps:

1. Run: no-subject-save2folder.jl:

In this script, we collected all emails that had an empty subject line. For each one, we created a folder and saved a txt file with all of the email information so we could check it manually (original-allprint.txt), and txt file with abbreviated information you the sender/receiver/cc/time/subject (original-short.txt).

2. Manually check errors

We manually opened and checked each original-allprint.txt file and original-short.txt,  and manually looked up the corresponding email in the original fauci email pdf (leopold-nih-foia-anthony-fauci-emails.pdf) to compare it with. We then saved the corrected email body text in a new file "fixed-body.txt" and manually saved the corrected sender/receiver/cc/date/subject in a file "fixed-short.txt". In some cases, this involved manually looking up an individual's email in the "names" vector of the original JSON (e.g., if Patricia Condrad was supposed to be in the cc list but this was not properly processed, we looked up the fact that she is the fifth person in the names array and saved this in the fixed-short.txt file).

3. Manual processing errors in emails with non-empty subjects

See long-subject-save2folders.jl 

In order to identify problems with emails that did not have an empty subject line, we first printed out all of the subject lines that did not start with "RE", "Re", "FW" or "Fw" (as these indicate a properly processed subject line). We then manually inspected each subject line, and flagged 25 subject lines for further inspection as they looked suspicious and possibly erroneous. Checking against the original pdf, we found that only 11 of these in fact contained errors, many of them minor. 

For the 11 non-empty subject emails with errors, we performed the same manual corrections as in step 2. 

4. Saving manual corrections to new JSON file

See save-manual-fixes.jl

We loaded in the original JSON file, and added 6 new entries to the "names" vector, corresponding to individuals that appeared in only one email but were not captured in the original automated processing due to parsing errors.

In this script, we loaded in information from the "fixed-body.txt" and "fixed-short.txt" files for each manually corrected file. This was saved into the final cleaned dataset:

fauci-email-data.json

5. Checking

See check-manual-fixes.jl.

As an additional sanity check, we read through the changed entries of the fauci-email-data.json file and output the new body and sender/receiver/cc/date/subject as txt files in the folders that contained our manual fixes. The fixed JSON file entries can then be checked against the original JSON entries, and the manual fixes we identified, to confirm that our manual fixes were properly incorporated.


Error analysis
-----------------

Not all of the emails we flagged and manually inspected needed fixes. 

Of the 57 emails with no subject lines, we have the following rough counts:

* 29 had issues with the subject being moved to the email body, but there were no errors in the sender/receiver/cc information
* 12 emails had no error (e.g., the original email did have an empty subject line)
* 16 emails required adjustments in the sender/receiver/cc information.

For the 25 possibly problematic emails with longer subject lines, we found:

* 5 examples where the CC was incorrectly listed as empty, the real CC was moved to the subject, and the subject was moved to the email body text
* 2 examples where part of the "to" field was moved to the subject
* 3 examples where part or all of the subject was moved to the email body text
* 15 examples that did not have any errors in the subject or sender/receiver/cc information

This means that of 25 emails flagged as possibly problematic, less than half were in fact erroneous, and there were only 7 cases of errors that affected information about sender and recipients (rather than just subject line parsing errors).

Overall, after manually inspecting all emails with empty subject lines and all emails with longer subject lines not starting with common patterns (FW, Fw, Re, RE), we found only 23 emails where an error was made with respect to sender/receiver/cc information. Given that there are 2761 emails in the dataset, this means an error rate of only 0.83% for these types of errors. All of them are fixed in the updated dataset.


Errors that remain
-------------------

There remain two common types of errors: OCR errors, and errors with timestamps (e.g., in a number of places we noticed the automated processing incorrectly processed times as AM instead of PM). We note that the exact time-of-day information will be particularly challenging to process regardless of the approach, as the emails were sent from many different time zones and were often shown in different formats.



