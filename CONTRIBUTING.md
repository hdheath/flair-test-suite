## How to contribute to FLAIR

We're very glad you want to help improve our code. Here's how you can contribute:

**Did you find a bug?** Please [file a ticket](https://github.com/BrooksLabUCSC/flair/issues). Be sure to tell us how you installed FLAIR, and how to reproduce the error.

**Do you have a feature request?** Please tell us about it in our [github forum](https://github.com/BrooksLabUCSC/flair/discussions). This is also the place to suggest improvements to our [documentation](https://flair.readthedocs.io/en/latest/) if anything is unclear.

**Do you want to fix something yourself?** To make improvements to the code directly please fork the repository and add your changes, then create a pull request (don't forget to describe what you did). You can do this even if you're just fixing typos in our documentation, any help is greatly appreciated.

### Logging guidelines

- Use `logging.getLogger(__name__)` in each module and avoid direct use of the root logger.
- The project formatter records the logger name, so manual prefixes like `[stage]` are unnecessary.
- Reserve `click.echo` for messages intended for the console; use the logger for run-summary files.

Please don't hesitate to reach out to us via the [forum](https://github.com/BrooksLabUCSC/flair/discussions).
Thanks!

The FLAIR team
