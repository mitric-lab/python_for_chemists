## How to interact with this website

This section gives an introduction on how to interact with the lecture notes.
These are organized into *chapters*.
Each chapter is a separate page.
Chapters are nested into a hierarchy of sub-chapters.
Typically, each chapter will be organized into a series of *headings* to subdivide a chapter.

### Navigation

There are several methods for navigating through the chapters of a book.

The **sidebar** on the left provides a list of all chapters.
Clicking on any of the chapter titles will load that page.

The sidebar may not automatically appear if the window is too narrow, particularly on mobile displays.
In that situation, the menu icon (three horizontal bars) at the top-left of the page can be pressed to open and close the sidebar.

The **arrow buttons** at the bottom of the page can be used to navigate to the previous or the next chapter.

The **left and right arrow keys** on the keyboard can be used to navigate to the previous or the next chapter.

### Top menu bar

The menu bar at the top of the page provides some icons for interacting with the notes.

| Icon | Description |
|------|-------------|
| <i class="fa fa-bars"></i> | Opens and closes the chapter listing sidebar. |
| <i class="fa fa-paint-brush"></i> | Opens a picker to choose a different color theme. |
| <i class="fa fa-search"></i> | Opens a search bar for searching within the book. |
| <i class="fa fa-print"></i> | Instructs the web browser to print the set of notes. |

Tapping the menu bar will scroll the page to the top.

### Search

The lecture notes have a built-in search system.
Pressing the search icon (<i class="fa fa-search"></i>) in the menu bar or pressing the `S` key on the keyboard will open an input box for entering search terms.
Typing any terms will show matching chapters and sections in real time.

Clicking any of the results will jump to that section.
The up and down arrow keys can be used to navigate the results, and enter will open the highlighted section.

After loading a search result, the matching search terms will be highlighted in the text.
Clicking a highlighted word or pressing the `Esc` key will remove the highlighting.

### Code blocks

Code blocks contain a copy icon <i class="fa fa-copy"></i>, that copies the code block into your local clipboard. 

Here's an example:

```python
print("Hello, World!")
```

We will often use the `assert` statement in code listings to show you
the value of a variable. Since the code blocks in this document are not
interactive  (you can not simply execute them in your browser), it is
not possible to print the value of the variable to the screen. 
Therefore, we ensure for you that all code blocks in these lecture notes 
run without errors and in this way we can represent the value of a
variable through the use of `assert`. 
The following code block for example shows that the variable `a` has
the value 2:
```python
a = 2
assert a == 2
```
If the condition would evaluate to `False`, the `assert` statement would
raise an `AssertionError`. 
