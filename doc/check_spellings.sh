#!/bin/sh

IGNORE_LIST="aspell_ignore_list.txt"

if ! command -v aspell &>/dev/null; then
    echo "Error: aspell is not installed. Please install it first."
    exit 1
fi

NEW_MISSPELLED_WORDS=$(find . -name "*.rst" | xargs cat | aspell list --lang=en_US --mode=sgml | sort | uniq)

UNSEEN_WORDS=$(echo "$NEW_MISSPELLED_WORDS" | grep -Fxvf "$IGNORE_LIST")

if [ -z "$UNSEEN_WORDS" ]; then
    echo "No new misspelled words found."
    echo "Spell check completed."
    exit 0
fi

echo "New misspelled words found:"
echo "$UNSEEN_WORDS"

for WORD in $UNSEEN_WORDS; do
    read -p "Add '$WORD' to the ignore list? (y/n): " RESPONSE
    if [[ "$RESPONSE" =~ ^[Yy]$ ]]; then
        echo "$WORD" >> "$IGNORE_LIST"
        sort -o "$IGNORE_LIST" "$IGNORE_LIST"  # Keep the list sorted
        echo "'$WORD' added to ignore list."
    fi
done

echo "Spell check completed."
