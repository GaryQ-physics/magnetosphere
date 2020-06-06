def urlretrieve(url, fname):
    """Python 2/3 urlretrieve compatability function.

    If Python 3, returns
    urllib.request.urlretrieve(url, fname)

    If Python 2, returns
    urllib.urlretrieve(url, fname)
    """

    import sys

    if sys.version_info[0] > 2:
        import urllib.request, urllib.error
        try:
            res = urllib.request.urlretrieve(url, fname)
            return res
        except urllib.error.URLError as e:
            print(e)
        except ValueError as e:
            print("'" + url + "' is not a valid URL")
    else:
        import urllib, urllib2, ssl
        try:
            context = ssl._create_unverified_context()
            res = urllib.urlretrieve(url, fname, context=context)
            return res
        except urllib2.URLError as e:
            print(e)
        except ValueError:
            print("'" + url + "' is not a valid URL")
